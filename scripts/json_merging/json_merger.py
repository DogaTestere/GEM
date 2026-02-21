#!/usr/bin/env python3

import json
import sqlite3
import argparse
import sys

from pathlib import Path
from datetime import datetime, timezone

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

def main(args):

    db_path = Path(args.kegg_db)
    json_path = Path(args.parsed_json)

    if not db_path.exists():
        raise FileNotFoundError(f"Database not found: {db_path}")

    if not json_path.exists():
        raise FileNotFoundError(f"Gene JSON not found: {json_path}")

    with json_path.open("r", encoding="utf-8") as f:
        gene_json = json.load(f)

    conn = connect_db(db_path)

    create_gene_table(conn)
    insert_gene_data(conn, gene_json)

    # ---- OPTIONAL GO ----
    if args.go_json:
        go_path = Path(args.go_json)

        if not go_path.exists():
            raise FileNotFoundError(f"GO JSON not found: {go_path}")

        with go_path.open("r", encoding="utf-8") as f:
            go_json = json.load(f)

        print("GO JSON detected — currently not processed.")
        # Future: insert_go_data(conn, go_json)

    conn.close()

    print("Gene metadata successfully added.")

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
    parser.add_argument("--parsed_json", required=True)
    parser.add_argument("--go_json", required=False)
    parser.add_argument("--kegg_db", required=True)
    args = parser.parse_args()

    main(args)

    write_versions()