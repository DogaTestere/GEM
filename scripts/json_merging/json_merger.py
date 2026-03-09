#!/usr/bin/env python3

import json
import sqlite3
import argparse
import sys
import shutil

from pathlib import Path
from datetime import datetime, timezone
from goatools.obo_parser import GODag

# --- GO ID for Transport Reaction Finding ---
# GO:0006810 <-- ID for 'Transport'
# GO:0005215 <-- ID for 'transporter activity'
# It climbs up the tree that belongs to these id's and tags the ones that uses it
# It also tries to find where it happens with cellular_compenent

GO_TRANSPORT_ROOTS = {
    "GO:0006810",  
    "GO:0005215",  
}

GO_CELLULAR_COMPONENT = "cellular_component"

class GOTransportClassifier:
    def __init__(self, go_obo_path):
        self.go_dag = GODag(go_obo_path)

        self.transport_terms = self._get_all_descendants(GO_TRANSPORT_ROOTS)

        self.compartment_map = {
            go_id: term.name
            for go_id, term in self.go_dag.items()
            if term.namespace == GO_CELLULAR_COMPONENT
        }

    def get_compartment(self, go_id):
        return self.compartment_map.get(go_id)
    
    def _get_all_descendants(self, root_ids):
        """
        Given a set of GO root IDs, return a set containing:
        - the roots
        - all descendant GO terms (recursive)
        """

        descendants = set()
        stack = list(root_ids)

        while stack:
            go_id = stack.pop()

            if go_id in descendants:
                continue

            descendants.add(go_id)

            term = self.go_dag.get(go_id)
            if not term:
                continue

            for child in term.children:
                stack.append(child.id)

        return descendants
    
# --- Class Ended ---
    
def connect_db(db_path):
    return sqlite3.connect(db_path)

def create_gene_table(conn):
    cur = conn.cursor()

    cur.execute("""
    CREATE TABLE IF NOT EXISTS genes (
        uniprot_id TEXT PRIMARY KEY,
        gene_name TEXT,
        protein_name TEXT,
        ec_number TEXT
    )
    """)

    conn.commit()

def create_go_table(conn):
    cur = conn.cursor()

    cur.execute("""
    CREATE TABLE IF NOT EXISTS go_terms (
        uniprot_id TEXT,
        go_id TEXT,
        go_name TEXT,
        qualifier TEXT,
        aspect TEXT,
        evidence TEXT,
        PRIMARY KEY (uniprot_id, go_id)
    )
    """)

    conn.commit()
    
    # Index for go terms
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_go_terms_uniprot
        ON go_terms(uniprot_id);
    """)
    conn.commit()
    
    cur.execute("""
    CREATE INDEX IF NOT EXISTS idx_gene_map_reaction
    ON gene_map (reaction_id);
    """)
    conn.commit()

def insert_gene_data(conn, gene_json):

    cur = conn.cursor()

    for entry in gene_json:

        uniprot = entry.get("UniProt_ID")
        gene = entry.get("Gene")
        protein = entry.get("Protein")
        ec = entry.get("EC_Number")

        if not uniprot:
            continue

        cur.execute("""
            INSERT INTO genes (uniprot_id, gene_name, protein_name, ec_number)
            VALUES (?, ?, ?, ?)
            ON CONFLICT(uniprot_id) DO UPDATE SET
                gene_name=excluded.gene_name,
                protein_name=excluded.protein_name,
                ec_number=excluded.ec_number
        """, (uniprot, gene, protein, ec))

    conn.commit()

def insert_go_data(conn, go_json):
    cur = conn.cursor()

    for uniprot_id, annotations in go_json.items():
        for entry in annotations:
            go_id = entry.get("goId")
            go_name = entry.get("goName")
            qualifier = entry.get("qualifier")
            aspect = entry.get("goAspect")
            evidence = entry.get("goEvidence")

            if not go_id:
                continue

            cur.execute("""
                INSERT INTO go_terms (uniprot_id, go_id, go_name, qualifier, aspect, evidence)
                VALUES (?, ?, ?, ?, ?, ?)
                ON CONFLICT(uniprot_id, go_id) DO NOTHING
            """, (uniprot_id, go_id, go_name, qualifier, aspect, evidence))

    conn.commit()
    
def load_gene_reaction_map(conn):
    cur = conn.cursor()
    cur.execute("SELECT uniprot_id, reaction_id FROM gene_map")

    mapping = {}
    for uniprot_id, reaction_id in cur.fetchall():
        mapping.setdefault(uniprot_id, []).append(reaction_id)

    return mapping

def annotate_transport_reactions(conn, go_json, classifier):
    cur = conn.cursor()

    # preload gene->reaction map
    gene_reaction_map = load_gene_reaction_map(conn)

    gene_compartments = {}
    transport_genes = set()

    # classify genes
    for uniprot_id, annotations in go_json.items():
        for entry in annotations:
            go_id = entry.get("goId")
            
            if not go_id:
                continue

            if go_id in classifier.transport_terms:
                transport_genes.add(uniprot_id)

            compartment = classifier.compartment_map.get(go_id)
            if compartment:
                existing = gene_compartments.get(uniprot_id)
                if existing is None or len(compartment) > len(existing):
                    gene_compartments[uniprot_id] = compartment

    # batch update reactions
    updates = []

    for gene in transport_genes:

        compartment = gene_compartments.get(gene)

        for reaction_id in gene_reaction_map.get(gene, []):
            updates.append((compartment, reaction_id))

    cur.executemany("""
        UPDATE reactions
        SET is_transport=1,
            compartment=COALESCE(?, compartment)
        WHERE reaction_id=?
    """, updates)

    conn.commit()

def main(
    input_db_path,
    output_db_path,
    gene_json_path,
    go_json_path=None,
    go_basic_path=None
):

    input_db = Path(input_db_path)
    output_db = Path(output_db_path)
    gene_json_path = Path(gene_json_path)

    if not input_db.exists():
        raise FileNotFoundError(input_db)

    if not gene_json_path.exists():
        raise FileNotFoundError(gene_json_path)

    shutil.copy(input_db, output_db)

    conn = connect_db(output_db)

    # ---- Insert gene metadata ----
    with gene_json_path.open("r", encoding="utf-8") as f:
        gene_json = json.load(f)

    create_gene_table(conn)
    insert_gene_data(conn, gene_json)

    # ---- GO processing ----
    if go_json_path and go_basic_path:

        go_json_path = Path(go_json_path)
        go_basic_path = Path(go_basic_path)

        if not go_json_path.exists():
            raise FileNotFoundError(go_json_path)

        if not go_basic_path.exists():
            raise FileNotFoundError(go_basic_path)

        print("Loading GO DAG...")

        classifier = GOTransportClassifier(str(go_basic_path))

        with go_json_path.open("r", encoding="utf-8") as f:
            go_json = json.load(f)

        create_go_table(conn)
        insert_go_data(conn, go_json)

        annotate_transport_reactions(conn, go_json, classifier)

        print("GO transport annotation complete.")

    conn.close()

    print("Database successfully created:", output_db)
    
def write_versions():
    versions = {
        "json_merging": {
            "python": sys.version.split()[0],
            "timestamp": datetime.now(timezone.utc).isoformat()
        }
    }

    with open("versions.yml", "w") as f:
        json.dump(versions, f, indent=2)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--input_db", required=True)
    parser.add_argument("--output_db", required=True)

    parser.add_argument("--parsed_json", required=True)

    parser.add_argument("--go_json")
    parser.add_argument("--go_basic")

    args = parser.parse_args()

    main(
        input_db_path=args.input_db,
        output_db_path=args.output_db,
        gene_json_path=args.parsed_json,
        go_json_path=args.go_json,
        go_basic_path=args.go_basic
    )

    write_versions()