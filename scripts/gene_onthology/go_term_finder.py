#!/usr/bin/env python3

import sys
import json
import argparse

from datetime import datetime, timezone
from bioservices import QuickGO

def extract_uniprot_ids(parsed_json):
    """
    Extract unique UniProt IDs from parsed GBK JSON.
    """
    return {
        entry["UniProt_ID"]
        for entry in parsed_json
        if entry.get("UniProt_ID")
    }


def fetch_go_terms(uniprot_ids):
    """
    Uses QuickGO from bioservices to get go names and go aspects from each uniprot id
    """
    go = QuickGO()
    go_data = {}

    for uid in uniprot_ids:
        try:
            response = go.Annotation(
                geneProductId=uid,        
                includeFields="goName"    
            )

            go_data[uid] = response.get("results", [])

        except Exception as e:
            sys.stderr.write(f"GO error for {uid}: {e}\n")
            go_data[uid] = []

    return go_data

def write_json(data, output_path):
    """
        Saves the term finding results as a json. Structure is this: \n
        {
            "P00561": [
                {
                    "id": "UniProtKB:P00561!612814768",
                    "geneProductId": "UniProtKB:P00561",
                    "qualifier": "acts_upstream_of_or_within",
                    "goId": "GO:0009089",
                    "goName": "L-lysine biosynthetic process via diaminopimelate",
                    "goEvidence": "IDA",
                    "goAspect": "biological_process",
                    "evidenceCode": "ECO:0000314",
                    "reference": "PMID:11368768",
                    "withFrom": null,
                    "taxonId": 83333,
                    "taxonName": null,
                    "assignedBy": "EcoCyc",
                    "extensions": null,
                    "targetSets": null,
                    "symbol": "thrA",
                    "date": "20070802",
                    "synonyms": null,
                    "name": null
                },
    """
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)

def write_versions():
    versions = {
        "json_merging": {
            "python": sys.version.split()[0],
            "timestamp": datetime.now(timezone.utc).isoformat()
        }
    }

    with open("versions.yml", "w") as f:
        json.dump(versions, f, indent=2)

def main():
    parser = argparse.ArgumentParser(
        description="Fetch GO terms for UniProt IDs using QuickGO"
    )
    parser.add_argument("--input_json", required=True)
    parser.add_argument("--go_json", required=True)

    args = parser.parse_args()

    with open(args.input_json) as f:
        parsed_gbk = json.load(f)

    uniprot_ids = extract_uniprot_ids(parsed_gbk)

    if not uniprot_ids:
        sys.stderr.write("No UniProt IDs found in input JSON\n")

    go_data = fetch_go_terms(uniprot_ids)

    write_json(go_data, args.go_json)
    write_versions()


if __name__ == "__main__":
    main()
