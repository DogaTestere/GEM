#!/usr/bin/env python3

import sqlite3
import os
import time
import requests
import sys
import json
import argparse
import re

from pathlib import Path
from datetime import datetime, timezone

BASE_URL = "https://rest.kegg.jp"

# --- SQL Database Creation ---

def open_kegg_db(db_path):
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON;")
    return conn

def create_tables(conn: sqlite3.Connection) -> None:
    """Create all required tables and indexes."""
    schema = """
    CREATE TABLE IF NOT EXISTS ec_numbers (
        uniprot_id TEXT,
        ec_number TEXT,
        sysname TEXT,
        orthology TEXT,
        all_reac TEXT,
        reac TEXT,
        substrate TEXT,
        product TEXT,
        manual_check INTEGER
    );

    CREATE TABLE IF NOT EXISTS reactions (
        reaction_id TEXT PRIMARY KEY,
        name TEXT,
        definition TEXT,
        equation TEXT,
        stoichiometry TEXT,
        pathway TEXT,
        orthology TEXT,
        gene_name TEXT,
        manual_flag INTEGER
    );

    CREATE TABLE IF NOT EXISTS compounds (
        compound_id TEXT PRIMARY KEY,
        name TEXT,
        formula TEXT,
        reaction TEXT,       -- comma-separated list of R-numbers
        enzyme TEXT,         -- comma-separated list of EC numbers
        manual_flag INTEGER
    );

    CREATE TABLE IF NOT EXISTS manual_checks (
        entity_type TEXT,
        entity_id   TEXT,
        timestamp   DATETIME DEFAULT CURRENT_TIMESTAMP,
        note        TEXT
    );

    CREATE TABLE IF NOT EXISTS reaction_stoichiometry (
        reaction_id TEXT,
        compound_id TEXT,
        coeff REAL,
        PRIMARY KEY (reaction_id, compound_id),
        FOREIGN KEY (reaction_id) REFERENCES reactions(reaction_id),
        FOREIGN KEY (compound_id) REFERENCES compounds(compound_id)
    );

    CREATE TABLE IF NOT EXISTS pathway_reaction (
        pathway_name  TEXT,
        reaction_id TEXT,
        PRIMARY KEY (pathway_name, reaction_id),
        FOREIGN KEY (reaction_id) REFERENCES reactions(reaction_id)
    );

    CREATE INDEX IF NOT EXISTS idx_ec_uniprot      ON ec_numbers(uniprot_id);
    CREATE INDEX IF NOT EXISTS idx_reaction_id     ON reactions(reaction_id);
    CREATE INDEX IF NOT EXISTS idx_compound_id     ON compounds(compound_id);
    CREATE INDEX IF NOT EXISTS idx_stoich_reaction ON reaction_stoichiometry(reaction_id);
    CREATE INDEX IF NOT EXISTS idx_stoich_compound ON reaction_stoichiometry(compound_id);
    CREATE INDEX IF NOT EXISTS idx_pathway_reaction ON pathway_reaction(pathway_name);
    """
    conn.executescript(schema)


def insert_data(conn, uniprot_id, gene_name, result_payload, manual_log):
    """
    Persist fully resolved KEGG entries to the database.
    Handles multiple reactions, compounds, and stoichiometries per EC number.
    """
    cursor = conn.cursor()

    ec_info = result_payload.get("ec_info")
    ec_number = result_payload.get("ec_number")
    reaction_infos = result_payload.get("reaction_infos", [])
    compound_info_dict = result_payload.get("compound_info", {})
    stio_infos = result_payload.get("stio_infos", {})

    # 1. Insert EC Info
    if ec_info:
        cursor.execute("""
            INSERT OR IGNORE INTO ec_numbers
            (uniprot_id, ec_number, sysname, orthology, all_reac, reac, substrate, product, manual_check)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            uniprot_id,
            ec_number,                                        
            ec_info.get("sysname"),
            ec_info.get("ec_orthology"),
            ",".join(ec_info.get("all_reactions", [])),
            ec_info.get("reaction"),
            ",".join(ec_info.get("substrates", [])),
            ",".join(ec_info.get("products",   [])),
            int(ec_info.get("manual_check", True)),
        ))
        
        if ec_info.get("manual_check"):
            _log_manual_check(
                conn, cursor,
                entity_type="ec_number",
                entity_id=ec_number,
                manual_log=manual_log,
                note=f"Empty ALL_REAC for UniProt {uniprot_id}"
            )

    # 2. Insert Compounds (Do this before stoichiometry to satisfy foreign keys)
    if compound_info_dict:
        for cmp_id, cmp_info in compound_info_dict.items():
            needs_review = not cmp_id.startswith("C")       
            cursor.execute("""
                INSERT OR IGNORE INTO compounds
                (compound_id, name, formula, reaction, enzyme, manual_flag)
                VALUES (?, ?, ?, ?, ?, ?)
            """, (
                cmp_id,
                cmp_info.get("compound_name"),
                cmp_info.get("formula"),
                ",".join(cmp_info.get("reaction_involved", [])),
                ",".join(cmp_info.get("enzymes_involved",  [])),
                int(needs_review),                            
            ))

            if needs_review:
                _log_manual_check(
                    conn, cursor,
                    entity_type="compound",
                    entity_id=cmp_id,
                    manual_log=manual_log,
                    note=f"Non-standard compound ID (not C-prefixed); UniProt {uniprot_id}"
                )

    # 3. Insert Reactions & Stoichiometry & Pathways
    for rn_info in reaction_infos:
        if not rn_info:
            continue
            
        rn_id = rn_info.get("rn_id")
        if not rn_id:
            continue
            
        stoichiometry = stio_infos.get(rn_id, {})

        cursor.execute("""
            INSERT OR IGNORE INTO reactions
            (reaction_id, name, definition, equation, stoichiometry,
             pathway, orthology, manual_flag, gene_name)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            rn_id,                                               
            rn_info.get("reaction_name"),                        
            rn_info.get("readable_eq"),
            rn_info.get("equation_raw"),
            json.dumps(stoichiometry) if stoichiometry else None,
            rn_info.get("pathway"),
            rn_info.get("rn_orthology"),
            int(rn_info.get("manual_check", True)),
            gene_name
        ))

        # Insert Pathway linking
        pathway_name = rn_info.get("pathway")
        if pathway_name:
            cursor.execute("""
                INSERT OR IGNORE INTO pathway_reaction (pathway_name, reaction_id)
                VALUES (?, ?)
            """, (pathway_name, rn_id))

        # Insert Stoichiometry linking
        if stoichiometry:
            for cmp_id, coeff in stoichiometry.items():
                is_nonstandard = not cmp_id.startswith("C")
                # Ensure compound exists just in case it wasn't in compound_info_dict
                cursor.execute("""
                    INSERT OR IGNORE INTO compounds (compound_id, manual_flag)
                    VALUES (?, ?)
                """, (cmp_id, int(is_nonstandard)))

                cursor.execute("""
                    INSERT OR IGNORE INTO reaction_stoichiometry
                    (reaction_id, compound_id, coeff)
                    VALUES (?, ?, ?)
                """, (rn_id, cmp_id, coeff))              

        # Manual check for the reaction itself
        if rn_info.get("manual_check"):
            _log_manual_check(
                conn, cursor,
                entity_type="reaction",
                entity_id=rn_id,
                manual_log=manual_log,
                note=f"Reaction contains non-standard compound IDs or missing equation; gene {gene_name}, UniProt {uniprot_id}"
            )

    conn.commit()


def _log_manual_check(conn, cursor, entity_type, entity_id, manual_log, note=""):
    """
    Insert a row into manual_checks and (optionally) append to the JSON log list.
    Duplicate (entity_type, entity_id) pairs are silently ignored in the DB.
    """
    ts = datetime.now(timezone.utc).isoformat()

    cursor.execute("""
        INSERT OR IGNORE INTO manual_checks (entity_type, entity_id, timestamp, note)
        VALUES (?, ?, ?, ?)
    """, (entity_type, entity_id, ts, note))

    if manual_log is not None:
        manual_log.append({
            "entity_type": entity_type,
            "entity_id":   entity_id,
            "timestamp":   ts,
            "note":        note,
        })

# ---------------------------------------------------------------------------
# KEGG API Fetching Helpers
# ---------------------------------------------------------------------------

def kegg_request(endpoint, rate_limit=0.35, retries=5, backoff=1.5, last_time=[0]):
    url = f"{BASE_URL}{endpoint}"

    for attempt in range(1, retries + 1):
        elapsed = time.time() - last_time[0]
        if elapsed < rate_limit:
            time.sleep(rate_limit - elapsed)

        try:
            print(f"      → GET {endpoint} (attempt {attempt})", file=sys.stderr)
            r = requests.get(url, timeout=(5, 20))
            r.raise_for_status()
            last_time[0] = time.time()
            return r.text

        except requests.RequestException as e:
            print(
                f"      ! ERROR {endpoint} attempt {attempt}: {e}",
                file=sys.stderr
            )
            if attempt == retries:
                return None
            time.sleep(backoff ** attempt)


def parse_page_entries(page, section):
    """
    Returns the content lines belonging to *section* from a KEGG flat-file page.
    Section headers are left-aligned uppercase words; continuation lines are
    indented by 12 spaces.
    """
    if not page:
        return []

    lines  = []
    active = False
    indent = " " * 12

    for line in page.splitlines():
        if line and not line.startswith(" "):
            if line.startswith(section):
                active  = True
                content = line.replace(section, "", 1).strip()
                if content:
                    lines.append(content)
            else:
                if active:
                    break
                active = False
        elif active:
            if line.startswith(indent):
                lines.append(line[12:].rstrip())

    return lines


# ---------------------------------------------------------------------------
# KEGG Information Parsers
# ---------------------------------------------------------------------------

def get_ec_info(pageResults):
    """
    Parses /get/{EC_Number}.
    Which is like this:
        ENTRY       EC 3.5.99.6                 Enzyme
        NAME        glucosamine-6-phosphate deaminase;
                    glucosaminephosphate isomerase (ambiguous);
                    glucosamine-6-phosphate isomerase (ambiguous);
                    phosphoglucosaminisomerase (ambiguous);
                    glucosamine phosphate deaminase;
                    aminodeoxyglucosephosphate isomerase (ambiguous);
                    phosphoglucosamine isomerase (ambiguous);
                    2-amino-2-deoxy-D-glucose-6-phosphate aminohydrolase (ketol isomerizing)
        CLASS       Hydrolases;
                    Acting on carbon-nitrogen bonds, other than peptide bonds;
                    In other compounds
        SYSNAME     2-amino-2-deoxy-alpha-D-glucose-6-phosphate aminohydrolase (ketol isomerizing)
        REACTION    alpha-D-glucosamine 6-phosphate + H2O = D-fructose 6-phosphate + NH3 [RN:R00765]
        ALL_REAC    R00765
        SUBSTRATE   alpha-D-glucosamine 6-phosphate;
                    H2O [CPD:C00001]
        PRODUCT     D-fructose 6-phosphate [CPD:C00085];
                    NH3 [CPD:C00014]
        ...
        PATHWAY     ec00520  Amino sugar and nucleotide sugar metabolism
                    ec01100  Metabolic pathways
        ORTHOLOGY   K02564  glucosamine-6-phosphate deaminase

    Returns
    -------
    dict with keys:
        reaction_id   : first R-number from REACTION line (e.g. "R00765")
        reaction      : full REACTION line string
        substrates    : [CPD-id, ...]
        products      : [CPD-id, ...]
        sysname       : systematic enzyme name string
        all_reactions : [R-number, ...]
        ec_orthology  : KO id (e.g. "K02564")
        manual_check  : True when there are empty returns of ALL_REAC

    Note: Some EC returns return a reaction with a macro molecule(ex. 2.1.1.198), thus needs to be marked for manual review. This is also used to find obsolete enzymes(ex. 4.99.1.1)
    """
    if not pageResults:
        return None

    sysname_lines  = parse_page_entries(pageResults, "SYSNAME")
    reaction_lines = parse_page_entries(pageResults, "REACTION")
    all_reac_lines = parse_page_entries(pageResults, "ALL_REAC")
    substrate_lines= parse_page_entries(pageResults, "SUBSTRATE")
    product_lines  = parse_page_entries(pageResults, "PRODUCT")
    orthology_lines= parse_page_entries(pageResults, "ORTHOLOGY")

    sysname = " ".join(sysname_lines).strip() if sysname_lines else None

    reaction_line = None
    reaction_id   = None
    for line in reaction_lines:
        match = re.search(r"\[RN:(R\d+)\]", line)
        if match:
            reaction_id   = match.group(1)
            reaction_line = line
            break

    all_reactions = []
    for line in all_reac_lines:
        all_reactions.extend(re.findall(r"R\d+", line))

    def extract_compounds(lines):
        compounds = []
        for line in lines:
            compounds.extend(re.findall(r"\[CPD:(C\d+)\]", line))
        return compounds

    substrates = extract_compounds(substrate_lines)
    products   = extract_compounds(product_lines)

    ko_id = None
    if orthology_lines:
        match = re.search(r"(K\d+)", orthology_lines[0])
        if match:
            ko_id = match.group(1)

    manual_check = len(all_reactions) == 0

    return {
        "reaction_id":   reaction_id,    
        "reaction":      reaction_line,
        "substrates":    substrates,
        "products":      products,
        "sysname":       sysname,
        "all_reactions": all_reactions,
        "ec_orthology":  ko_id,
        "manual_check":  manual_check
    }


def get_reaction_info(pageResults):
    """
    Parses /get/{R-number}.
    Which is like this:
        ENTRY       R00765                      Reaction
        NAME        D-glucosamine-6-phosphate aminohydrolase (ketol isomerizing);
                    D-glucosamine-6-phosphate ketol-isomerase(deaminating)
        DEFINITION  D-Glucosamine 6-phosphate + H2O <=> D-Fructose 6-phosphate + Ammonia
        EQUATION    C00352 + C00001 <=> C00085 + C00014
        RCLASS      RC00163  C00085_C00352
        ENZYME      3.5.99.6
        PATHWAY     rn00520  Amino sugar and nucleotide sugar metabolism
        ...
        ORTHOLOGY   K02564  glucosamine-6-phosphate deaminase [EC:3.5.99.6]
        DBLINKS     RHEA: 12175

    Returns
    -------
    dict with keys:
        reaction_name : human-readable NAME string (If NAME isn't avaliable it fetches COMMENT instead, and if both are unavaliable it uses the DEFINITION as a name)
        rn_id         : R-number itself (Which can be parsed from ENTRY)
        readable_eq   : DEFINITION string
        equation_raw  : raw EQUATION string (compound IDs)
        substrates    : ordered list of (compound_id, coeff) tuples
        products      : ordered list of (compound_id, coeff) tuples
        pathway       : human readable part of PATHWAY string
        rn_orthology  : KO id
        manual_check  : True when reaction_id and/or equation is empty
        
    Note : Some reactions contain macromolecules or non-C compounds(ex. G10481(n+2)) these are marked with manual_check
    """
    if not pageResults:
        return None

    entry_lines      = parse_page_entries(pageResults, "ENTRY")
    name_lines       = parse_page_entries(pageResults, "NAME")
    comment_lines    = parse_page_entries(pageResults, "COMMENT")
    definition_lines = parse_page_entries(pageResults, "DEFINITION")
    equation_lines   = parse_page_entries(pageResults, "EQUATION")
    pathway_lines    = parse_page_entries(pageResults, "PATHWAY")
    orthology_lines  = parse_page_entries(pageResults, "ORTHOLOGY")

    rn_id = None
    if entry_lines:
        match = re.match(r"(R\d+)", entry_lines[0].strip())
        if match:
            rn_id = match.group(1)

    readable_eq = " ".join(definition_lines) if definition_lines else None
    equation    = " ".join(equation_lines)   if equation_lines   else None

    if name_lines:
        reaction_name = name_lines[0].rstrip(";")
    elif comment_lines:
        reaction_name = comment_lines[0]
    else:
        reaction_name = readable_eq 

    substrates   = []
    products     = []
    manual_check = not bool(rn_id and equation)

    if equation:
        if "<=>" in equation:
            left, right = equation.split("<=>", 1)
        elif "=>" in equation:
            left, right = equation.split("=>", 1)
        elif "<=" in equation:
            right, left = equation.split("<=", 1)
        else:
            left = right = None

        def parse_side(side):
            nonlocal manual_check 
            ordered = []
            if not side:
                return ordered
            tokens = re.split(r'\s+\+\s+', side.strip())
            for token in tokens:
                token = token.strip()
                match = re.match(r'^(\d+(?:\.\d+)?)?\s*([A-Z]\d+(?:\([^)]*\))?)', token)
                if match:
                    coeff_str = match.group(1)
                    raw_id    = match.group(2)
                    coeff     = float(coeff_str) if coeff_str else 1.0
                    clean_id  = re.sub(r'\([^)]*\)', '', raw_id)
                    ordered.append((clean_id, coeff))
                    
                    if not clean_id.startswith("C"):
                        manual_check = True
            return ordered

        substrates = parse_side(left)
        products   = parse_side(right)

    pathway_name = None
    if pathway_lines:
        match = re.search(r"rn\d+\s+(.*)", pathway_lines[0])
        if match:
            pathway_name = match.group(1).strip()
        else:
            pathway_name = pathway_lines[0].strip()


    ko_id = None
    if orthology_lines:
        match = re.search(r"(K\d+)", orthology_lines[0])
        if match:
            ko_id = match.group(1)

    return {
        "rn_id":         rn_id,           
        "reaction_name": reaction_name,   
        "readable_eq":   readable_eq,
        "equation_raw":  equation,
        "substrates":    substrates,
        "products":      products,
        "pathway":       pathway_name,
        "rn_orthology":  ko_id,
        "manual_check":  manual_check,
    }


def get_compound_info(pageResults):
    """
    Parses /get/{compound_id}.
    Which looks like this:
        ENTRY       C00352                      Compound
        NAME        D-Glucosamine 6-phosphate;
                    D-Glucosamine phosphate
        FORMULA     C6H14NO8P
        ...
        REACTION    R00765 R00768 R01961 R01964 R01965 R02058 R02059 R02060 
                    R02631 R12220 R12594
        ...
        ENZYME      2.3.1.4         2.6.1.16        2.7.1.1         2.7.1.8         
                    2.7.1.147       2.7.1.-         3.5.1.25        3.5.99.6        
                    5.1.3.42        5.4.2.10
    Returns
    -------
    dict with keys:
        compound_name     : NAME string
        formula           : molecular formula (may be None)
        reaction_involved : [R-number, ...]
        enzymes_involved  : [EC-number, ...]
        manual_check      : True if the ENTRY part contains non-C compound_id
    """
    if not pageResults:
        return None
    
    entry_lines    = parse_page_entries(pageResults, "ENTRY")
    name_lines     = parse_page_entries(pageResults, "NAME")
    formula_lines  = parse_page_entries(pageResults, "FORMULA")
    reaction_lines = parse_page_entries(pageResults, "REACTION")
    enzyme_lines   = parse_page_entries(pageResults, "ENZYME")

    compound_name = name_lines[0].rstrip(";") if name_lines else None
    formula       = formula_lines[0]          if formula_lines else None

    reaction_list = []
    for line in reaction_lines:
        reaction_list.extend(re.findall(r"R\d+", line))

    enzyme_list = []
    for line in enzyme_lines:
        enzyme_list.extend(re.findall(r"\d+\.\d+\.\d+\.\d+|\d+\.\d+\.\d+\.-", line))
        
    manual_check = False
    if entry_lines:
        # Match standard KEGG IDs (1 letter followed by digits, e.g., C00352, G10481)
        match = re.match(r"^([A-Za-z]\d+)", entry_lines[0].strip())
        if match:
            compound_id = match.group(1)
            # Flag if the ID does not start with 'C'
            if not compound_id.startswith("C"):
                manual_check = True
        else:
            # Flag if an ID couldn't be cleanly parsed
            manual_check = True
    else:
        # Flag if the ENTRY field is completely missing
        manual_check = True

    return {
        "compound_name":     compound_name,
        "formula":           formula,
        "reaction_involved": reaction_list,
        "enzymes_involved":  enzyme_list,
        "manual_check":      manual_check
    }

# ---------------------------------------------------------------------------
# Functions for Unfinished / Incomplete EC Numbers
# ---------------------------------------------------------------------------

def get_unfinished_ecs(pageResults):
    """
    Parses a CDS page for entries with incomplete EC numbers (e.g. 3.5.99.-) from /get/{kegg_id}.
    Which looks like this:
    ENTRY       JW3105            CDS       T00068
    SYMBOL      agaS
    NAME        (GenBank) tagatose-6-phosphate ketose/aldose isomerase
    ORTHOLOGY   K02082  D-galactosamine 6-phosphate deaminase/isomerase [EC:3.5.99.-]
    ORGANISM    ecj  Escherichia coli K-12 W3110
    PATHWAY     ecj00052  Galactose metabolism
                ecj01100  Metabolic pathways
    ...

    Returns
    -------
    dict with keys:
        unfinished_ko : K-number string
        gene_id       : gene symbol string
        cds_name      : gene name string 
        org_id        : organism string 
    """
    if not pageResults:
        return None

    symbol_lines   = parse_page_entries(pageResults, "SYMBOL")
    name_lines     = parse_page_entries(pageResults, "NAME")
    orthology_lines= parse_page_entries(pageResults, "ORTHOLOGY")
    organism_lines = parse_page_entries(pageResults, "ORGANISM")

    ko_id = None
    if orthology_lines:
        match = re.search(r"(K\d+)", orthology_lines[0])
        if match:
            ko_id = match.group(1)

    return {
        "unfinished_ko": ko_id,
        "gene_id":       symbol_lines[0]   if symbol_lines   else None,
        "cds_name":      name_lines[0]     if name_lines     else None,
        "org_id":        organism_lines[0] if organism_lines else None,
    }


def get_rn_id(pageResults):
    """
    Parses a /link/rn/{ko_id} page to return the reaction connected to it.
    Which looks like this:
    ko:K12524	rn:R00480
    ko:K12524	rn:R01773
    ko:K12524	rn:R01775

    Returns
    -------
    list: [Rxxxx, Rxxx, ...]
    """
    if not pageResults:
        return []

    result = []
    
    lines = pageResults.splitlines() if isinstance(pageResults, str) else pageResults

    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        match = re.search(r"rn:(R\d+)", line)
        if match:
            rn_id = match.group(1)
            
            if rn_id not in result:
                result.append(rn_id)

    return result


# ---------------------------------------------------------------------------
# Stoichiometry Builder
# ---------------------------------------------------------------------------

def generate_stio(reaction_info):
    """
    This uses the output from get_reaction_info(). It skips writing stio info for anything marked with manual_check. (It is just a precaution to filling wrong stio info before manual check fixes)
    
    It uses the substrates and products values with their coeffiecents from the output.
    
    Returns
    ---------
    dict   { compound_id : coeff }
    """
    stoichiometry = {}

    if not reaction_info or reaction_info.get("manual_check"):
        return stoichiometry

    for cid, coeff in reaction_info.get("substrates", []):
        stoichiometry[cid] = stoichiometry.get(cid, 0.0) - coeff

    for cid, coeff in reaction_info.get("products", []):
        stoichiometry[cid] = stoichiometry.get(cid, 0.0) + coeff

    return stoichiometry


# ---------------------------------------------------------------------------
# Module-level caches — persist for the lifetime of a single process run
# ---------------------------------------------------------------------------
reaction_cache = {}
compound_cache = {}
ec_cache       = {}
ko_cache       = {}


# ---------------------------------------------------------------------------
# Main Orchestrator
# ---------------------------------------------------------------------------

def run_kegg(conn, ec_number=None, kegg_id=None, gene_name=None):
    ec_info        = None
    reaction_infos = []
    rn_ids         = []
    link_info      = [] 
    
    incomplete_ec = False

    if ec_number in ec_cache:
        ec_info = ec_cache[ec_number]
        rn_id = ec_info.get("reaction_id") if ec_info else None
        rn_ids = [rn_id] if rn_id else []
    else:
        ec_page = kegg_request(f"/get/{ec_number}")
        
        if ec_page is None: 
            ec_info = None
            incomplete_ec = True
        else:
            ec_info = get_ec_info(ec_page)
            ec_cache[ec_number] = ec_info
            if ec_info:
                rn_id = ec_info.get("reaction_id")
                rn_ids = [rn_id] if rn_id else []
            
    if ec_info is None and incomplete_ec == True:
        cds_page = kegg_request(f"/get/{kegg_id}")
        
        if cds_page:
            cds_info = get_unfinished_ecs(cds_page)
            ko_id = cds_info.get("unfinished_ko") if cds_info else None
            
            if ko_id:
                if ko_id in ko_cache:
                    link_info = ko_cache[ko_id]
                else:
                    ko_rn_link_page = kegg_request(f"/link/rn/{ko_id}")
                    if ko_rn_link_page:
                        link_info = get_rn_id(ko_rn_link_page)
                        ko_cache[ko_id] = link_info
            
            if link_info:
                rn_ids = link_info 
    
    for rn_id in rn_ids:
        if rn_id in reaction_cache:
            reaction_infos.append(reaction_cache[rn_id])
        else:
            reaction_page = kegg_request(f"/get/{rn_id}")
            if reaction_page:
                rn_info = get_reaction_info(reaction_page)
                reaction_cache[rn_id] = rn_info
                reaction_infos.append(rn_info)
            
    compound_list = []
    for rn_info in reaction_infos:
        if rn_info:
            compound_list += [cid for cid, _ in rn_info.get("substrates", [])]
            compound_list += [cid for cid, _ in rn_info.get("products", [])]
    
    for comp_id in compound_list:
        if comp_id in compound_cache:
            compound_info = compound_cache[comp_id]
        else:
            compound_page = kegg_request(f"/get/{comp_id}")
            if compound_page:
                compound_info = get_compound_info(compound_page)
                compound_cache[comp_id] = compound_info
    
        
    stio_tables = {}
    for rn_info in reaction_infos:
        if rn_info and rn_info.get("rn_id"):
            stio_tables[rn_info["rn_id"]] = generate_stio(rn_info)
            
    return {
        "ec_info"       : ec_info,
        "ec_number"     : ec_number,
        "reaction_infos": reaction_infos,
        "rn_ids"        : rn_ids,
        "all_linked_rns": link_info, 
        "compound_info" : {cid: compound_cache[cid] for cid in compound_list if cid in compound_cache},
        "stio_infos"    : stio_tables,
        "gene_name"     : gene_name
    }   



# ---------------------------------------------------------------------------
# Entry Point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch KEGG data from JSON mapping")
    parser.add_argument("--db",            required=True)
    parser.add_argument("--parsed_json",   required=True)
    parser.add_argument("--kegg_map_json", required=True)
    parser.add_argument("--manual_json",   required=True)
    args = parser.parse_args()

    conn = open_kegg_db(args.db)
    create_tables(conn)

    with open(args.parsed_json) as f:
        parsed_entries = json.load(f)
    with open(args.kegg_map_json) as f:
        kegg_map = json.load(f)

    manual_log = []

    for entry in parsed_entries:
        uni_id    = entry.get("UniProt_ID")
        ec_number = entry.get("EC_Number")
        kegg_id   = kegg_map.get(uni_id)
        gene_name = entry.get("Gene")

        if not ec_number and not kegg_id:
            continue

        result = run_kegg(
            conn,
            ec_number=ec_number,
            kegg_id=kegg_id,
            gene_name=gene_name,
        )

        if result:
            insert_data(
                conn,
                uniprot_id=uni_id,
                gene_name=gene_name,
                result_payload=result,
                manual_log=manual_log
            )

    # ------------------------------------------------------------------
    # Write JSON manual-check log
    # ------------------------------------------------------------------
    manual_json_path = Path(args.manual_json)
    with open(manual_json_path, "w") as f:
        json.dump(
            {
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "source_db":    args.db,
                "total_flags":  len(manual_log),
                "entries":      manual_log,
            },
            f,
            indent=2,
        )
    print(f"Manual-check log written → {manual_json_path}", file=sys.stderr)
    print("All entries processed successfully.", file=sys.stderr)


