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

def open_kegg_db(db_path):
    """
    Creates or opens a SQLite database for KEGG information.
    Refuses creation inside Nextflow work directory.
    """

    db_path = Path(db_path).expanduser()

    if not db_path.is_absolute():
        db_path = db_path.resolve()

    # Prevent DB inside Nextflow work/
    nxf_work = os.environ.get("NXF_WORK")
    if nxf_work:
        nxf_work = Path(nxf_work).resolve()
        try:
            db_path.relative_to(nxf_work)
            raise RuntimeError(
                f"Refusing to create persistent DB inside Nextflow work dir: {db_path}"
            )
        except ValueError:
            pass

    db_path.parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(db_path, timeout=30)

    cur = conn.cursor()
    try:
        cur.execute("PRAGMA journal_mode=WAL;")
        cur.execute("PRAGMA synchronous=NORMAL;")
    except sqlite3.DatabaseError:
        pass

    return conn


def create_tables(conn):
    cur = conn.cursor()

    # --- Reactions ---
    cur.execute("""
    CREATE TABLE IF NOT EXISTS reactions (
        reaction_id TEXT PRIMARY KEY,
        reaction_name TEXT,
        stoichiometry TEXT,
        lower_bound REAL,
        upper_bound REAL
    )
    """)

    # --- Compounds ---
    cur.execute("""
    CREATE TABLE IF NOT EXISTS compounds (
        compound_id TEXT PRIMARY KEY,
        name TEXT,
        formula TEXT,
        dblinks TEXT
    )
    """)

    # --- Reaction-Compound map ---
    cur.execute("""
    CREATE TABLE IF NOT EXISTS reaction_compounds (
        reaction_id TEXT,
        compound_id TEXT,
        stoichiometry REAL,
        PRIMARY KEY (reaction_id, compound_id)
    )
    """)

    # ---Gene map ---
    cur.execute("""
    CREATE TABLE IF NOT EXISTS gene_map (
        uniprot_id TEXT,
        kegg_id TEXT,
        reaction_id TEXT,
        PRIMARY KEY (uniprot_id, kegg_id, reaction_id)
    )
    """)

    # --- Completion tracking ---
    cur.execute("""
    CREATE TABLE IF NOT EXISTS kegg_status (
        kegg_id TEXT PRIMARY KEY,
        completed INTEGER
    )
    """)

    # ---------------- Indexes ----------------

    # Fast lookup of reactions by ID
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_reactions_id
        ON reactions (reaction_id)
    """)

    # Fast lookup of compounds by ID
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_compounds_id
        ON compounds (compound_id)
    """)

    # Fast lookup of reaction → compound mapping
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_reaction_compounds_reaction
        ON reaction_compounds (reaction_id)
    """)

    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_reaction_compounds_compound
        ON reaction_compounds (compound_id)
    """)

    # Fast lookup for gene-based queries
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_gene_map_uniprot
        ON gene_map (uniprot_id)
    """)

    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_gene_map_kegg
        ON gene_map (kegg_id)
    """)

    conn.commit()


def is_kegg_completed(conn, kegg_id):
    cur = conn.cursor()
    cur.execute(
        "SELECT completed FROM kegg_status WHERE kegg_id=?",
        (kegg_id,)
    )
    row = cur.fetchone()
    return row is not None and row[0] == 1

def mark_kegg_completed(conn, kegg_id):
    cur = conn.cursor()
    cur.execute("""
        INSERT OR REPLACE INTO kegg_status (kegg_id, completed)
        VALUES (?, 1)
    """, (kegg_id,))
    conn.commit()

def insert_reaction(conn, reaction_id, reaction_data):
    cur = conn.cursor()

    cur.execute("""
        INSERT OR IGNORE INTO reactions
        (reaction_id, reaction_name, lower_bound, upper_bound)
        VALUES (?, ?, ?, ?)
    """, (
        reaction_id,
        reaction_data["reaction_name"],
        reaction_data["lower_bound"],
        reaction_data["upper_bound"]
    ))

    conn.commit()

def insert_compound(conn, compound_id, compound_data):
    cur = conn.cursor()

    cur.execute("""
        INSERT OR IGNORE INTO compounds
        (compound_id, name, formula, dblinks)
        VALUES (?, ?, ?, ?)
    """, (
        compound_id,
        compound_data["name"],
        compound_data["formula"],
        json.dumps(compound_data["dblinks"])
    ))

    conn.commit()

def insert_reaction_compounds(conn, reaction_id, stoichiometry_dict):
    cur = conn.cursor()

    for compound_id, coeff in stoichiometry_dict.items():
        cur.execute("""
            INSERT OR REPLACE INTO reaction_compounds
            (reaction_id, compound_id, stoichiometry)
            VALUES (?, ?, ?)
        """, (reaction_id, compound_id, coeff))

    conn.commit()


def insert_gene_map(conn, uniprot_id, kegg_id, reaction_id):
    cur = conn.cursor()
    cur.execute("""
        INSERT OR IGNORE INTO gene_map
        VALUES (?, ?, ?)
    """, (uniprot_id, kegg_id, reaction_id))
    conn.commit()

# --- KEGG API FETCHİNG ---

def kegg_request(endpoint, rate_limit=0.35, retries=5, backoff=1.5, last_time=[0]):
    url = f"{BASE_URL}{endpoint}"

    for attempt in range(1, retries + 1):
        elapsed = time.time() - last_time[0]
        if elapsed < rate_limit:
            time.sleep(rate_limit - elapsed)

        try:
            print(f"      → GET {endpoint} (attempt {attempt})", file=sys.stderr)
            r = requests.get(url, timeout=(5, 20))  # CONNECT, READ
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
    Returns the whole page request with each section seperated
    """

    if not page:
        return []

    lines = []
    active = False
    # All KEGG entries align section names. Content starts 12 spaces in.
    indent = " " * 12

    for line in page.splitlines():
        # Section headers are always left-aligned, uppercase
        if line and not line.startswith(" "):
            if line.startswith(section):
                active = True
                # Section headers may have content on same line
                content = line.replace(section, "", 1).strip()
                if content:
                    lines.append(content)
            else:
                # New section started; if we were reading ours, stop
                if active:
                    break
                active = False
        elif active:
            # Continuation lines start with indent
            if line.startswith(indent):
                lines.append(line[12:].rstrip())

    return lines

def find_brite_numbers(pageResult):
    """
    Find BRITE part's numbers \n
    Similar to find_ko_numbers() but with (ec\d{5})

    Example page: \n
    BRITE       KEGG Orthology (KO) [BR:ecc00001] \n
    ... \n
    Enzymes [BR:ecc01000] \n
    ...
    """
    brites = parse_page_entries(pageResult, "BRITE")
    brite_ids = []
    brite_regex = re.compile(r'([a-z]{3,5}\d{5})')

    for line in brites:
        # Find the first ID on the line
        match = brite_regex.match(line)
        if match:
            brite_ids.append(match.group(1))

    return brite_ids


def find_full_brite(pageResult):
    """
    Finds the full BRITE section in the page (with the description) \n
        BRITE       KEGG Orthology (KO) [BR:ecc00001] \n
             09190 Not Included in Pathway or Brite \n
              09191 Unclassified: metabolism \n
               99980 Enzymes with EC numbers \n
                c3299 (ygbL) \n
            Enzymes [BR:ecc01000] \n
             4. Lyases \n
              4.1  Carbon-carbon lyases \n
               4.1.1  Carboxy-lyases \n
                4.1.1.104  3-dehydro-4-phosphotetronate decarboxylase \n
                 c3299 (ygbL) \n
    """
    brites = parse_page_entries(pageResult, "BRITE")
    return brites

def find_pathways(pageResults):
    """
    Searches the PATHWAY id and returns them. \n

    Example: \n
    PATHWAY     ecj00260  Glycine, serine and threonine metabolism \n
                    ecj01100  Metabolic pathways \n
                    ecj01110  Biosynthesis of secondary metabolites \n
                    ecj01230  Biosynthesis of amino acids \n
    
    Returns the whole ecj... included parts
    """

    pathways = parse_page_entries(pageResults, "PATHWAY")
    pathway_ids = []

    pathway_regex = re.compile(r'^([a-z]{3}\d{5})\b')

    for line in pathways:
        match = pathway_regex.match(line)
        if match:
            pathway_ids.append(match.group(1))

    return pathway_ids

def find_org_pathways(pageResults):
    """
    Searches the ecj... coded pathways to get the: \n
    GENE        JW0001  thrA; fused aspartokinase I and homoserine dehydrogenase I [KO:K12524] [EC:2.7.2.4 1.1.1.3] \n
    part.

    Then returns the thrA part and KO:..... part \n
    Returns list of dicts: \n
    [ \n
        {"gene": "thrA", "ko": "K12524"}, \n
        ... \n
    ]
    """
    genes = parse_page_entries(pageResults, "GENE")
    results = []

    ko_regex = re.compile(r'\[KO:(K\d{5})\]')
    gene_regex = re.compile(r'^\S+\s+([^;]+)')

    for line in genes:
        ko_match = ko_regex.search(line)
        gene_match = gene_regex.search(line)

        if ko_match and gene_match:
            gene_name = gene_match.group(1).strip()
            ko_id = ko_match.group(1)

            results.append({
                "gene": gene_name,
                "ko": ko_id
            })

    return results

def find_ko_rn_links(pageResults):
    """
    Reads the result of /link/rn/ko:Kxxxx and returns the rn:Rxxxx information \n
    Example: \n
        ko:K12524	rn:R00480 \n
        ko:K12524	rn:R01773 \n
        ko:K12524	rn:R01775 \n
    """
    if not pageResults:
        return []

    rn_ids = []
    rn_regex = re.compile(r'rn:(R\d{5})')

    for line in pageResults.splitlines():
        match = rn_regex.search(line)
        if match:
            rn_ids.append(match.group(1))

    rn_ids = list(dict.fromkeys(rn_ids))
    return rn_ids

def find_reactions(pageResults):
    """
    Reads get/rn:Rxxx page and finds \n
    NAME        ATP:L-aspartate 4-phosphotransferase \n
    DEFINITION  ATP + L-Aspartate <=> ADP + 4-Phospho-L-aspartate \n
    EQUATION    C00002 + C00049 <=> C00008 + C03082 \n

    It seperates the compund ids while preserving order. It also seperates the arrow indicator. \n
    <= : Right side is reactants, saved as (-1000,0) \n
    <=> : Left side is reactants, saved as (-1000,1000) \n
    => : Left side is reactants, saved as (0,1000)
    """

    name = " ".join(parse_page_entries(pageResults, "NAME"))
    equation = " ".join(parse_page_entries(pageResults, "EQUATION"))

    reactants = []
    products = []

    # Determine bounds
    if "<=>" in equation:
        lower_bound, upper_bound = -1000.0, 1000.0
        left, right = equation.split("<=>", 1)
    elif "=>" in equation:
        lower_bound, upper_bound = 0.0, 1000.0
        left, right = equation.split("=>", 1)
    elif "<=" in equation:
        lower_bound, upper_bound = -1000.0, 0.0
        right, left = equation.split("<=", 1)
    else:
        lower_bound, upper_bound = 0.0, 1000.0
        left = equation
        right = ""

    react_dict, react_list = parse_side(left, -1)
    prod_dict, prod_list = parse_side(right, +1)

    stoichiometry = {**react_dict, **prod_dict}

    # Preserve order + remove duplicates
    reactants = list(dict.fromkeys(reactants))
    products = list(dict.fromkeys(products))

    return {
        "reaction_name": name,
        "reactants": react_list,
        "products": prod_list,
        "stoichiometry": stoichiometry,
        "lower_bound": lower_bound,
        "upper_bound": upper_bound
    }

def parse_side(side, sign):
    compounds = {}
    ordered_ids = []

    parts = side.split("+")
    for part in parts:
        part = part.strip()

        match = re.match(r'(?:(\d+)\s*)?(C\d{5})', part)
        if match:
            coeff = int(match.group(1)) if match.group(1) else 1
            cid = match.group(2)

            compounds[cid] = sign * coeff
            ordered_ids.append(cid)

    return compounds, ordered_ids

def find_compunds(pageResults):
    """
    Reads get/Cxxxx page and finds \n
    NAME        ATP; \n
                Adenosine 5'-triphosphate \n
    FORMULA     C10H16N5O13P3 \n
    ... \n
    DBLINKS     CAS: 56-65-5 \n
            PubChem: 3304 \n
            ChEBI: 15422 \n
            KNApSAcK: C00001491 \n
            PDB-CCD: ATP \n
            NIKKAJI: J10.680A \n
    
    It saves them as a list, in case of multiple names, only the first name is saved.
    """
    names = parse_page_entries(pageResults, "NAME")
    formula_lines = parse_page_entries(pageResults, "FORMULA")
    dblinks_lines = parse_page_entries(pageResults, "DBLINKS")

    compound_name = None
    if names:
        compound_name = names[0].split(";")[0].strip()

    formula = formula_lines[0] if formula_lines else None

    dblinks = {}
    for line in dblinks_lines:
        if ":" in line:
            key, val = line.split(":", 1)
            dblinks[key.strip()] = val.strip()

    return {
        "name": compound_name,
        "formula": formula,
        "dblinks": dblinks
    }

def fetch_kegg_results(kegg_id):
    """
    Calls other functions and organises the results.
    """
    ko_cache = {}
    reaction_cache = {}
    compound_cache = {}

    results = {
        "BRITE_ids" : [],
        "BRITE" : [],
        "Pathway_ids" : [],
        "Pathway_KOs" : [],
        "Reaction_ids" : [],
        "Reaction_results" : [],
        "Compound_results" : [],
    }

    print("[ENTRY-MAIN] Fetching main page", file=sys.stderr)
    page = kegg_request(f"/get/{kegg_id}")
    if not page:
        print("[ENTRY-ERROR] Failed on page retrival", file=sys.stderr)
        return results
    
    print("[ENTRY-BRITE] Finding BRITE information from main page", file=sys.stderr)
    results["BRITE"] = find_full_brite(page)
    results["BRITE_ids"] = find_brite_numbers(page)

    print("[ENTRY-PATH] Finding Pathway information from main page", file=sys.stderr)
    pathways = find_pathways(page)
    results["Pathway_ids"] = pathways

    # --- Pathway KO numbers ---
    # {"gene": "thrA", "ko": "K12524"} for these

    all_KOs = []

    for pathway_id in pathways:
        print(f"[PATH-MAIN] Fetching pathway page, id:{pathway_id}", file=sys.stderr)
        pathway_page = kegg_request(f"/get/{pathway_id}")
        if not pathway_page:
            print(f"[PATH-ERROR] Fetching pathway page, id:{pathway_id}", file=sys.stderr)

        print("[PATH-KO] Fetching KO ids", file=sys.stderr)
        org_KOs = find_org_pathways(pathway_page)
        for item in org_KOs:
            ko_id = item["ko"]
            if ko_id not in all_KOs:
                all_KOs.append(ko_id)
    
    results["Pathway_KOs"] = all_KOs

    # --- KO ---> RN_ID transformation ---
    # ko:K12524	rn:R00480 for these

    all_rns = []

    for ko_id in all_KOs:
        if ko_id in ko_cache:
            rn_ids = ko_cache[ko_id]
        else:
            print(f"[KO-MAIN] Fetching ko-rn page, id:{ko_id}", file=sys.stderr)
            ko_page = kegg_request(f"/link/rn/ko:{ko_id}")

            print("[KO-LINK] Fetching reaction ids", file=sys.stderr)
            rn_ids = find_ko_rn_links(ko_page)
            ko_cache[ko_id] = rn_ids
        
        for rn_id in rn_ids:
            if rn_id not in all_rns:
                all_rns.append(rn_id)
    
    # --- Reaction info ---
    # DEFINITION  ATP + L-Aspartate <=> ADP + 4-Phospho-L-aspartate for these
    all_comp = set()

    for rn_id in all_rns:
        if rn_id in reaction_cache:
            reaction_info = reaction_cache[rn_id]
        else:
            print(f"[RN-MAIN] Fetching reaction page, id:{rn_id}")
            rn_page = kegg_request(f"/get/{rn_id}")
            if not rn_page:
                continue

            reaction_info = find_reactions(rn_page)
            reaction_cache[rn_id] = reaction_info

        reaction_info_with_id = {
            "reaction_id": rn_id,
            **reaction_info
        }

        results["Reaction_results"].append(reaction_info_with_id)

        for comp in reaction_info["reactants"] + reaction_info["products"]:
            all_comp.add(comp)

    # --- Compound info ---
    # FORMULA     C10H16N5O13P3 for these

    for comp_id in all_comp:
        if comp_id in compound_cache:
            comp_info = compound_cache[comp_id]
        else:
            print(f"[COMP-MAIN] Fetching compund page, id:{comp_id}")
            comp_page = kegg_request(f"/get/{comp_id}")
            if not comp_page:
                continue

            comp_info = find_compunds(comp_page)
            compound_cache[comp_id] = comp_info

        results["Compound_results"].append({
            "compound_id": comp_id,
            **comp_info
        })

    return results
    
def write_versions():
    versions = {
        "json_merging": {
            "python": sys.version.split()[0],
            "timestamp": datetime.now(timezone.utc).isoformat()
        }
    }

    with open("versions.yml", "w") as f:
        json.dump(versions, f, indent=2)

def run_db_kegg_fetch(mapping_json, db_path):
    """
    Reads a mapping json file like the following: \n
    "A0A0H2VA12": "ecc:c3299", \n
    "A0A0H2VA68": "ecc:c3297", \n
    "A0A0H2VDN9": "ecc:c5321", \n
    Right side is UniProt id, left side is kegg_ids \n

    Then checks the db, if there is no db or info is missing it runs fetch_kegg_results() to complete it.
    Returns just the db
    """
    with open(mapping_json) as f:
        mapping = json.load(f)

    # ---------------- Connect DB ----------------
    with open(mapping_json) as f:
        mapping = json.load(f)

    conn = open_kegg_db(db_path)
    create_tables(conn)

    # ---------------- Process each UniProt ----------------
    for uniprot_id, kegg_id in mapping.items():

        print(f"[MAIN] Processing {uniprot_id} -> {kegg_id}",
              file=sys.stderr)

        if is_kegg_completed(conn, kegg_id):
            print("   → Already completed. Skipping.",
                  file=sys.stderr)
            continue

        # Fetch from KEGG
        kegg_data = fetch_kegg_results(kegg_id)

        # ---------------- Store reactions ----------------
        for reaction in kegg_data["Reaction_results"]:

            reaction_id = reaction["reaction_id"]

            insert_reaction(conn, reaction_id, reaction)

            insert_reaction_compounds(conn, reaction_id, reaction["stoichiometry"])

            insert_gene_map(conn, uniprot_id, kegg_id, reaction_id)

        # ---------------- Store compounds ----------------
        for compound in kegg_data["Compound_results"]:
            insert_compound(
                conn,
                compound["compound_id"],
                compound
            )

        mark_kegg_completed(conn, kegg_id)

    conn.close()

    return db_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping_json", required=True)
    parser.add_argument("--db", required=True)
    args = parser.parse_args()

    run_db_kegg_fetch(
        mapping_json=args.mapping_json,
        db_path=args.db,
    )

    write_versions()