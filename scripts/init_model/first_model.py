import argparse
import json
import sys
import re
import sqlite3
import logging
from collections import defaultdict

from cobra import Model, Reaction, Metabolite
from cobra.io import write_sbml_model, validate_sbml_model

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

def SId_check(raw_id) -> str:
    """Make a string SBML compliant (SId type)."""
    if raw_id is None:
        return None
    raw_id = str(raw_id)
    clean = re.sub(r'[^a-zA-Z0-9_]', '_', raw_id)
    if not re.match(r'[a-zA-Z_]', clean):
        clean = '_' + clean
    return clean

def compartment_code(compartment_name: str) -> str:
    mapping = {
        'cytosol': 'c', 'cytoplasm': 'c', 'extracellular': 'e', 'periplasm': 'p',
        'mitochondrion': 'm', 'nucleus': 'n', 'endoplasmic reticulum': 'er',
        'golgi apparatus': 'g', 'peroxisome': 'x', 'vacuole': 'v'
    }
    name_lower = compartment_name.strip().lower()
    return mapping.get(name_lower, name_lower[:3])

# ------------------------------------------------------------------------------
# Database loading & Parsing
# ------------------------------------------------------------------------------

def load_database(db_path: str):
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    # Compounds (Bileşik isimleri ve formülleri için)
    compounds = {}
    try:
        cursor.execute("SELECT compound_id, name, formula FROM compounds")
        for row in cursor.fetchall():
            compounds[row['compound_id']] = dict(row)
    except sqlite3.OperationalError:
        logger.warning("Compounds tablosu bulunamadı, bileşik isimleri ID olarak bırakılacak.")

    # Reactions ve Stoichiometry (Aynı tablodan okunuyor)
    cursor.execute("""
        SELECT reaction_id, name, definition, stoichiometry, lower_bound, upper_bound, gene_name
        FROM reactions
    """)
    
    reactions = {}
    stoich_dict = {} # reaction_id -> {compound_id: katsayı}
    gene_rows = []

    for row in cursor.fetchall():
        rid = row['reaction_id']
        
        # Transport tespiti (İsminde veya tanımında geçiyorsa)
        text_for_search = f"{row['name']} {row['definition']}".lower()
        is_transport = any(word in text_for_search for word in ['transport', 'exchange', 'uptake', 'efflux'])
        
        reactions[rid] = {
            'name': row['name'] or rid,
            'lower_bound': row['lower_bound'] if row['lower_bound'] is not None else -1000.0,
            'upper_bound': row['upper_bound'] if row['upper_bound'] is not None else 1000.0,
            'is_transport': is_transport,
            'compartment': 'cytosol' # Tabloda olmadığı için hücre içi kabul ediyoruz
        }

        # JSON formatındaki Stokiyometriyi okuma
        stoich_str = row['stoichiometry']
        if stoich_str:
            try:
                parsed_stoich = json.loads(stoich_str)
                # JSON sözlüğünü Python sözlüğüne çevir { 'C00019': -1.0, ... }
                stoich_dict[rid] = {k: float(v) for k, v in parsed_stoich.items()}
            except json.JSONDecodeError:
                logger.error(f"Reaksiyon {rid} için stoichiometry JSON parse edilemedi!")
                stoich_dict[rid] = {}
        else:
            stoich_dict[rid] = {}

        # GPR (Gen Kuralları)
        if row['gene_name']:
            gene_rows.append((rid, row['gene_name']))

    conn.close()
    return compounds, reactions, stoich_dict, gene_rows

# ------------------------------------------------------------------------------
# Model Kurulumu (Metabolit ve Reaksiyonlar)
# ------------------------------------------------------------------------------

def sanitize_formula(formula):
    if formula is None: return None
    return re.sub(r'[\(\)\-nR]', '', str(formula))

def create_metabolites(compounds, reactions, stoich_dict):
    compound_compartments = defaultdict(set)
    reaction_compound_comp = {}

    for rid, metabolites in stoich_dict.items():
        if rid not in reactions: continue
        comp_name = reactions[rid]['compartment']
        comp_code = compartment_code(comp_name)
        
        for cid in metabolites.keys():
            compound_compartments[cid].add(comp_code)
            reaction_compound_comp[(rid, cid)] = comp_code

    metabolite_objs = {}
    for cid, comps in compound_compartments.items():
        base_info = compounds.get(cid, {'name': cid, 'formula': None})
        base_name = base_info.get('name') or cid
        formula = base_info.get('formula')
        
        for comp in comps:
            clean_id = SId_check(cid)
            met_id = f"{clean_id}_{comp}"
            met = Metabolite(
                id=met_id,
                name=f"{base_name} [{comp}]",
                formula=sanitize_formula(formula),
                compartment=comp
            )
            metabolite_objs[(cid, comp)] = met

    final_stoich = defaultdict(dict)
    for rid, metabolites in stoich_dict.items():
        for cid, coeff in metabolites.items():
            comp = reaction_compound_comp.get((rid, cid))
            met = metabolite_objs.get((cid, comp))
            if met:
                final_stoich[rid][met] = coeff

    return metabolite_objs, final_stoich

def build_gene_rules(gene_rows):
    gene_map = defaultdict(set)
    for rid, gene_name in gene_rows:
        if not gene_name: continue
        # Genler bazen virgül veya noktalı virgül ile ayrılmış olabilir
        genes = re.split(r'[;,]', gene_name)
        for g in genes:
            clean_gene = SId_check(g.strip())
            if clean_gene:
                gene_map[SId_check(rid)].add(clean_gene)

    rules = {}
    for rid, genes in gene_map.items():
        if len(genes) == 1:
            rules[rid] = f"({next(iter(genes))})"
        else:
            rules[rid] = "(" + " or ".join(sorted(genes)) + ")"
    return rules

def is_consumption(metabolite, coeff, lower, upper):
    if lower >= 0: return coeff < 0
    if upper <= 0: return coeff > 0
    return False

def create_reactions(reactions, stoich, gene_rules):
    reaction_objects = []

    for rid, info in reactions.items():
        clean_rid = SId_check(rid)
        name = info['name'] or clean_rid
        lower = info['lower_bound']
        upper = info['upper_bound']

        metabolites = stoich.get(rid, {})
        if not metabolites:
            continue

        is_simple_boundary = len(metabolites) == 1

        if is_simple_boundary:
            (met, coeff), = metabolites.items()
            comp = met.compartment
            reversible = (lower < 0 and upper > 0)

            if comp == 'e':
                rxn_type = 'exchange'
            else:
                if reversible:
                    rxn_type = 'sink'
                else:
                    if is_consumption(met, coeff, lower, upper):
                        rxn_type = 'demand'
                    else:
                        rxn_type = 'sink'
                        if lower >= 0: lower = -1000.0
                        if upper <= 0: upper = 1000.0

            rxn = Reaction(id=clean_rid, name=f"{rxn_type}_{name}", lower_bound=lower, upper_bound=upper)
            rxn.add_metabolites({met: coeff})
        else:
            rxn = Reaction(id=clean_rid, name=name, lower_bound=lower, upper_bound=upper)
            rxn.add_metabolites(metabolites)

        if clean_rid in gene_rules:
            rxn.gene_reaction_rule = gene_rules[clean_rid]

        reaction_objects.append(rxn)

    return reaction_objects

def build_model(db_path: str, model_name: str, objective: str = None, direction: str = "max"):
    logger.info("Loading data from database...")
    compounds, reactions, stoich_dict, gene_rows = load_database(db_path)

    logger.info("Creating metabolites...")
    metabolite_objs, stoich = create_metabolites(compounds, reactions, stoich_dict)

    logger.info("Building gene rules...")
    gene_rules = build_gene_rules(gene_rows)

    logger.info(f"Creating reactions... (Total found: {len(reactions)})")
    reaction_objs = create_reactions(reactions, stoich, gene_rules)

    model = Model(model_name)
    model.add_reactions(reaction_objs)

    if objective:
        clean_obj = SId_check(objective)
        if clean_obj in [r.id for r in model.reactions]:
            model.objective = clean_obj
            model.objective.direction = direction
            logger.info(f"Objective successfully set to {clean_obj}")
        else:
            if model.reactions:
                dummy_rxn = next(iter(model.reactions))
                model.objective = dummy_rxn
                model.objective.direction = direction
                logger.warning(f"Requested objective '{objective}' not found. Using '{dummy_rxn.id}' as dummy.")
            else:
                logger.error("Model has no reactions! Cannot set objective.")

    model.repair()
    return model


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build a COBRApy model from a metabolic database.")
    parser.add_argument("--model_objective", required=True, help="Reaction ID to set as objective")
    parser.add_argument("--model_direction", default="max", choices=["max", "min"], help="Objective direction")
    parser.add_argument("--balanced_db", required=True, help="Path to SQLite database")
    parser.add_argument("--output_model", required=True, help="Output SBML file name")
    args = parser.parse_args()

    model = build_model(
        db_path=args.balanced_db,
        model_name="first_model",
        objective=args.model_objective,
        direction=args.model_direction
    )

    logger.info(f"Writing model to {args.output_model}...")
    write_sbml_model(model, args.output_model)

    logger.info("Validating SBML...")
    _, report = validate_sbml_model(args.output_model)
    if any(report.values()):
        logger.warning("SBML validation warnings/errors:")
        for key, errors in report.items():
            if errors:
                logger.warning(f"{key}: {errors}")

    logger.info("Done.")