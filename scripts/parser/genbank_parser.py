#!/usr/bin/env python3

import argparse
import json
import sys
import tempfile
import re
import os
from datetime import datetime, timezone
from Bio import SeqIO, BiopythonWarning
import warnings

# Suppress BiopythonParserWarnings
warnings.simplefilter("ignore", BiopythonWarning)

def sanitize_genbank(in_path):
    """
    Fix non-standard LOCUS lines so Biopython can parse them.
    """
    fixed_lines = []
    with open(in_path, "r") as f:
        for line in f:
            if line.startswith("LOCUS"):
                # extract locus name (first non-space token after LOCUS)
                m_locus = re.match(r"LOCUS\s+(\S+)", line)
                locus = m_locus.group(1) if m_locus else "UNKNOWN"
                # try to find a bp length (digits before 'bp'), otherwise first integer in the line
                m_len = re.search(r"(\d+)\s*bp", line)
                if not m_len:
                    m_len = re.search(r"(\d+)", line)
                length = m_len.group(1) if m_len else "0"
                line = f"LOCUS       {locus:<16} {length} bp    DNA     linear   BCT 01-JAN-2000\n"
            fixed_lines.append(line)

    tmp = tempfile.NamedTemporaryFile(delete=False, mode="w")
    tmp.writelines(fixed_lines)
    tmp.close()
    return tmp.name

def parse_genbank(gb_file):
    """
    Extract CDS entries from a GenBank file.
    Returns two lists:
    1. EC-first: CDS with EC numbers
    2. UniProt-first: CDS with EC numbers and UniProt IDs
    """
    clean_file = sanitize_genbank(gb_file)

    ec_first_data = []
    uniprot_first_data = []
    total_cds = 0

    for record in SeqIO.parse(clean_file, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue

            total_cds += 1
            qualifiers = feature.qualifiers

            gene = qualifiers.get("gene", [None])[0]
            protein = qualifiers.get("product", [None])[0]
            ec_number = qualifiers.get("EC_number", [None])[0]

            # Skip CDS without EC numbers for EC-first list
            if not ec_number:
                continue

            # Add to EC-first list
            ec_first_data.append({
                "Gene": gene,
                "Protein": protein,
                "EC_Number": ec_number
            })

            # Check for UniProt ID
            uniprot_id = None
            for uni_ids in qualifiers.get("inference", []):
                if "UniProtKB:" in uni_ids:
                    uniprot_id = uni_ids.split("UniProtKB:")[-1]
                    break

            # Add to UniProt-first list if UniProt ID exists
            if uniprot_id:
                uniprot_first_data.append({
                    "Gene": gene,
                    "Protein": protein,
                    "EC_Number": ec_number,
                    "UniProt_ID": uniprot_id
                })

    if os.path.exists(clean_file):
        os.remove(clean_file)

    sys.stderr.write(f"Total CDS checked: {total_cds}\n")
    sys.stderr.write(f"EC-first records: {len(ec_first_data)}\n")
    sys.stderr.write(f"UniProt-first records: {len(uniprot_first_data)}\n")

    return ec_first_data, uniprot_first_data

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
    parser = argparse.ArgumentParser(description="Parse GenBank file and extract CDS entries")
    parser.add_argument("--gb_file", required=True, help="Input GenBank file (.gbk or .gff)")
    parser.add_argument("--ec_out", default="ec_first.json", help="Output JSON for CDS with EC numbers")
    parser.add_argument("--uniprot_out", default="uniprot_first.json", help="Output JSON for CDS with EC + UniProt IDs")
    args = parser.parse_args()

    ec_data, uniprot_data = parse_genbank(args.gb_file)

    with open(args.ec_out, "w") as f:
        json.dump(ec_data, f, indent=2)
    with open(args.uniprot_out, "w") as f:
        json.dump(uniprot_data, f, indent=2)

    sys.stderr.write(f"Saved {len(ec_data)} records to {args.ec_out}\n")
    sys.stderr.write(f"Saved {len(uniprot_data)} records to {args.uniprot_out}\n")

    write_versions()

if __name__ == "__main__":
    main()
