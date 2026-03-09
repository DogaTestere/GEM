import argparse
import json
import sys
import re
import sqlite3

from datetime import datetime, timezone
from collections import defaultdict
from cobra import Model, Reaction, Metabolite
from cobra.io import write_sbml_model, validate_sbml_model

def SId_check(raw_id:str) -> str:
    """
    This function checks if the id's are SBML compliant. Main check is that ids cannot start with numbers.
    It checks compound, reaction and gene id's \n
    letter   ::=   ’a’..’z’,’A’..’Z’ \n
    digit    ::=   ’0’..’9’ \n
    idChar   ::=   letter | digit | ’_’ \n
    SId      ::=   ( letter | ’_’ ) idChar* \n
    If it is not compliant, it changes it to be compliant.
    """
    clean = re.sub(r'[^a-zA-Z0-9_]', '_', raw_id)
    if not re.match(r'[a-zA-Z_]', clean):
        clean = '_' + clean
    return clean

def database_loading(db_path):
    """
    It reads the database and extracts the following tables: \n
    CREATE TABLE compounds ( \n
        compound_id TEXT PRIMARY KEY,\n
        name TEXT,\n
        formula TEXT,\n
        dblinks TEXT\n
    )\n
    CREATE TABLE gene_map (\n
        uniprot_id TEXT,\n
        kegg_id TEXT,\n
        reaction_id TEXT,\n
        PRIMARY KEY (uniprot_id, kegg_id, reaction_id)\n
    )\n
    CREATE TABLE genes (\n
        uniprot_id TEXT PRIMARY KEY,\n
        gene_name TEXT,\n
        protein_name TEXT,\n
        ec_number TEXT\n
    )\n
    CREATE TABLE reaction_compounds (\n
        reaction_id TEXT,\n
        compound_id TEXT,\n
        stoichiometry REAL,\n
        PRIMARY KEY (reaction_id, compound_id)\n
    )\n
    CREATE TABLE reactions (\n
        reaction_id TEXT PRIMARY KEY,\n
        reaction_name TEXT,\n
        lower_bound REAL,\n
        upper_bound REAL\n
    )
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # compounds
    compound_info = {
        row[0]: {"name": row[1], "formula": row[2]}
        for row in cursor.execute(
            "SELECT compound_id, name, formula FROM compounds"
        )
    }

    # reactions
    reaction_info = {
        row[0]: {
            "reaction_name": row[1],
            "lower_bound": row[2],
            "upper_bound": row[3],
        }
        for row in cursor.execute(
            "SELECT reaction_id, reaction_name, lower_bound, upper_bound FROM reactions"
        )
    }

    # stoichiometry (single pass, flat)
    stoich_rows = cursor.execute("""
        SELECT reaction_id, compound_id, stoichiometry
        FROM reaction_compounds
    """).fetchall()

    # genes (single pass)
    gene_rows = cursor.execute("""
        SELECT gm.reaction_id, g.gene_name
        FROM gene_map gm
        JOIN genes g ON gm.uniprot_id = g.uniprot_id
    """).fetchall()

    conn.close()

    return compound_info, reaction_info, stoich_rows, gene_rows

def create_metabolites(compound_info, compartment="c"):
    """
    Fills the Metabolite class with the values gathered from database: \n
    Metabolite(\n
        id='compound_id',\n
        formula='formula',\n
        name='name',\n
        compartment="c"\n
    )\n
    """
    metabolites = {}

    for raw_id, info in compound_info.items():

        clean = SId_check(raw_id)

        metabolites[raw_id] = Metabolite(
            id=f"{clean}[{compartment}]",
            name=info["name"] or clean,
            formula=info["formula"],
            compartment=compartment,
        )

    return metabolites
    
def build_stio(stoich_rows, metabolites):
    """
    Produces the following dict to later use in reaction creation:\n
    {\n
        Metabolite_object : coefficient\n
    }\n
    """

    stoich = defaultdict(dict)

    for reaction_id, compound_id, coeff in stoich_rows:

        stoich[SId_check(reaction_id)][metabolites[compound_id]] = coeff

    return stoich

def build_gene_info(gene_rows):
    """
    Produces the following string for later use in reaction creation: \n
    (geneA)\n
    (geneA or geneB)\n
    """
    gene_map = defaultdict(set)

    for reaction_id, gene_name in gene_rows:

        gene_map[SId_check(reaction_id)].add(SId_check(gene_name))

    rules = {}

    for rxn_id, genes in gene_map.items():

        if len(genes) == 1:
            rules[rxn_id] = f"({next(iter(genes))})"
        else:
            rules[rxn_id] = "(" + " or ".join(sorted(genes)) + ")"

    return rules

def create_reactions(reaction_info, stoich, gene_rules):
    """
    Fills the Reaction class with the values gathered from database: \n
    Reaction( \n
        id='reaction_id', \n
        name='reaction_name', \n
        subsystem= '', \n
        lower_bound='lower_bound', \n
        upper_bound='upper_bound' \n
    )

    reaction.gene_reaction_rule = genes \n
    reaction.add_metabolites(metabolites_dict) \n
    """
    reactions = []

    append = reactions.append 

    for raw_id, info in reaction_info.items():

        rxn_id = SId_check(raw_id)

        rxn = Reaction(
            id=rxn_id,
            name=info["reaction_name"] or rxn_id,
            lower_bound=info["lower_bound"] or -1000,
            upper_bound=info["upper_bound"] or 1000,
        )

        if rxn_id in stoich:
            rxn.add_metabolites(stoich[rxn_id])

        if rxn_id in gene_rules:
            rxn.gene_reaction_rule = gene_rules[rxn_id]

        append(rxn)

    return reactions

def build_model(db_path, model_name, objective=None, direction="max"):
    """
    This runs metabolite and reaction creation and setting of objective alongside direction
    If the model passes SMBL test then it is send to another scripts for balancing as xml model
    """

    compound_info, reaction_info, stoich_rows, gene_rows = database_loading(db_path)

    model = Model(model_name)

    metabolites = create_metabolites(compound_info)

    stoich = build_stio(stoich_rows, metabolites)

    gene_rules = build_gene_info(gene_rows)

    reactions = create_reactions(
        reaction_info,
        stoich,
        gene_rules
    )

    model.add_reactions(reactions)

    if objective:
        model.objective = SId_check(objective)
        model.objective.direction = direction

    model.repair()
    model.solver.update()

    return model

def write_versions():
    versions = {
        "first_model": {
            "python": sys.version.split()[0],
            "timestamp": datetime.now(timezone.utc).isoformat()
        }
    }

    with open("versions.yml", "w") as f:
        json.dump(versions, f, indent=2)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--model_objective", required=True)
    parser.add_argument("--model_direction", default="max")
    parser.add_argument("--balanced_db", required=True)
    parser.add_argument("--output_model", required=True)

    args = parser.parse_args()

    compound_info, reaction_compound_info, reaction_info, gene_map = database_loading(
        args.balanced_db
    )

    model = build_model(
        compound_info=compound_info,
        reaction_compound_info=reaction_compound_info,
        reaction_info=reaction_info,
        gene_map=gene_map,
        model_name="first_model",
        objective=args.model_objective,
        direction=args.model_direction
    )

    write_sbml_model(model, args.output_model)

    _, report = validate_sbml_model(args.output_model)

    if any(report.values()):
        raise RuntimeError(report)

    write_versions()