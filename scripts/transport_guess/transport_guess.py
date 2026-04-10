#!/usr/bin/env python3

import json
import sqlite3
import argparse
import shutil
from pathlib import Path
from goatools.obo_parser import GODag

# --- GO ID for Transport Reaction Finding ---
# GO:0006810 <-- ID for 'Transport'
# GO:0005215 <-- ID for 'transporter activity'
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

def create_transport_table(conn):
    """Creates the new isolated table and adds performance indices."""
    cur = conn.cursor()
    
    # Create the table
    cur.execute("""
    CREATE TABLE IF NOT EXISTS transport_annotations (
        uniprot_id TEXT,
        go_id TEXT,
        go_name TEXT,
        go_aspect TEXT,
        transport_location TEXT,
        transport_name TEXT
    )
    """)
    
    # Add indices mirroring the style of idx_ec_uniprot to speed up downstream Nextflow queries
    cur.execute("""
    CREATE INDEX IF NOT EXISTS idx_transport_uniprot 
    ON transport_annotations(uniprot_id)
    """)
    
    cur.execute("""
    CREATE INDEX IF NOT EXISTS idx_transport_goid 
    ON transport_annotations(go_id)
    """)
    
    conn.commit()

def process_transport_annotations(conn, go_json, classifier):
    cur = conn.cursor()

    # 1. Build a mapping of gene_name -> reaction name from the existing reactions table
    # This allows us to fill the "Full name of the transport" column
    cur.execute("SELECT gene_name, name FROM reactions WHERE gene_name IS NOT NULL")
    
    gene_to_reaction_name = {}
    for gene_name, rxn_name in cur.fetchall():
        # If a gene maps to multiple reactions, we combine their names
        if gene_name in gene_to_reaction_name:
            gene_to_reaction_name[gene_name].add(rxn_name)
        else:
            gene_to_reaction_name[gene_name] = {rxn_name}

    # 2. Find the best known location (compartment) for each gene
    gene_compartments = {}
    for uniprot_id, annotations in go_json.items():
        for entry in annotations:
            go_id = entry.get("goId")
            if not go_id:
                continue
            
            compartment = classifier.get_compartment(go_id)
            if compartment:
                existing = gene_compartments.get(uniprot_id)
                # Keep the longest string as a heuristic for the most descriptive compartment
                if existing is None or len(compartment) > len(existing):
                    gene_compartments[uniprot_id] = compartment

    # 3. Extract the transport records
    transport_records = []
    
    for uniprot_id, annotations in go_json.items():
        for entry in annotations:
            go_id = entry.get("goId")
            if not go_id:
                continue

            # Check if this specific GO ID is a transport term
            if go_id in classifier.transport_terms:
                go_name = entry.get("goName")
                go_aspect = entry.get("goAspect")
                symbol = entry.get("symbol") # e.g., "HI_0370" or "iscA"

                location = gene_compartments.get(uniprot_id, "Unknown")
                
                # Retrieve the full reaction name based on the gene symbol
                rxn_names = gene_to_reaction_name.get(symbol, {"Unknown"})
                transport_name = " | ".join(filter(None, rxn_names))

                transport_records.append((
                    uniprot_id, 
                    go_id, 
                    go_name, 
                    go_aspect, 
                    location, 
                    transport_name
                ))

    # 4. Insert data into our new table
    cur.executemany("""
        INSERT INTO transport_annotations (
            uniprot_id, go_id, go_name, go_aspect, transport_location, transport_name
        )
        VALUES (?, ?, ?, ?, ?, ?)
    """, transport_records)

    conn.commit()


def main(input_db_path, output_db_path, go_json_path, go_basic_path):
    input_db = Path(input_db_path)
    output_db = Path(output_db_path)
    go_json_path = Path(go_json_path)
    go_basic_path = Path(go_basic_path)

    if not input_db.exists():
        raise FileNotFoundError(f"Input DB missing: {input_db}")
    if not go_json_path.exists():
        raise FileNotFoundError(f"GO JSON missing: {go_json_path}")
    if not go_basic_path.exists():
        raise FileNotFoundError(f"GO OBO missing: {go_basic_path}")

    # Copy the database so we don't mutate the original inputs in the pipeline
    shutil.copy(input_db, output_db)
    conn = connect_db(output_db)

    print("Loading GO DAG...")
    classifier = GOTransportClassifier(str(go_basic_path))

    with go_json_path.open("r", encoding="utf-8") as f:
        go_json = json.load(f)

    print("Creating transport table and indices...")
    create_transport_table(conn)
    
    print("Processing and inserting transport annotations...")
    process_transport_annotations(conn, go_json, classifier)

    conn.close()
    print(f"Database successfully created with new transport_annotations table: {output_db}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find transport reactions and build an isolated annotations table.")
    
    parser.add_argument("--input_db", required=True)
    parser.add_argument("--output_db", required=True)
    parser.add_argument("--go_json", required=True)
    parser.add_argument("--go_obo", required=True)

    args = parser.parse_args()

    main(
        input_db_path=args.input_db,
        output_db_path=args.output_db,
        go_json_path=args.go_json,
        go_basic_path=args.go_obo
    )