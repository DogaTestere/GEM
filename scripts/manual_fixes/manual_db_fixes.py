import json
import sqlite3
import argparse
import shutil
import sys
import re

def parse_equation(equation: str):
    """Parses a KEGG-style reaction equation into substrates and products."""
    if not equation:
        return [], []
        
    if "<=>" in equation:
        left, right = equation.split("<=>", 1)
    elif "=>" in equation:
        left, right = equation.split("=>", 1)
    elif "<=" in equation:
        right, left = equation.split("<=", 1)
    else:
        left = right = None

    def parse_side(side):
        ordered = []
        if not side:
            return ordered
        tokens = re.split(r'\s+\+\s+', side.strip())
        for token in tokens:
            token = token.strip()
            if not token:
                continue
            match = re.match(r'^(\d+(?:\.\d+)?)?\s*([A-Z]\d+(?:\([^)]*\))?)', token)
            if match:
                coeff_str = match.group(1)
                raw_id    = match.group(2)
                coeff     = float(coeff_str) if coeff_str else 1.0
                clean_id  = re.sub(r'\([^)]*\)', '', raw_id)
                ordered.append((clean_id, coeff))
        return ordered

    substrates = parse_side(left)
    products   = parse_side(right)
    return substrates, products

def read_manual_fix(json_filepath: str, in_db_filepath: str, out_db_filepath: str) -> str:
    """
    Reads the manual fix json file and records them back to the database
    Fixes added should be in this format:
   {
      "entity_type": "compound",
      "entity_id": "G10609",
      "timestamp": "2026-03-09T19:42:08.689584+00:00",
      "note": "Non-standard compound ID (not C-prefixed); UniProt P77293",
      "fix" : {
          "formula" : "formula_here"
          },
    },
    {
      "entity_type": "reaction",
      "entity_id": "R12808",
      "timestamp": "2026-03-09T19:42:08.689728+00:00",
      "note": "Reaction contains non-standard compound IDs or missing equation; gene yfdH, UniProt P77293",
      "fix" : {
          "equation" : "equation_here"
          }
    },
    {
      "entity_type": "ec_number",
      "entity_id": "5.6.2.2",
      "timestamp": "2026-03-09T19:34:12.134266+00:00",
      "note": "Empty ALL_REAC for UniProt P0AES6",
      "fix" : {
          "all_reac" : "reaction_name",
          "missing_compounds" : {
              "substrates" : [
                  {
                      "name" : "compund_name",
                      "id" : "compund_id"
                  }
              ],
              "products" : [
                  {
                      "name" : "compound_name",
                      "id" : "compound_id"
                  }
              ]
         }
    },


    ! Currently ec_number fix doesn't work
    """ 
    try:
        shutil.copy2(in_db_filepath, out_db_filepath)
    except IOError as e:
        print(f"Error copying database: {e}")
        sys.exit(1)
        
    try:
        with open(json_filepath, 'r', encoding='utf-8') as f:
            fixes = json.load(f)
            
        # Failsafe: Standardize 'fixes' to always be a list of dictionaries
        if isinstance(fixes, dict):
            if "entity_type" in fixes:
                fixes = [fixes]
            else:
                extracted_list = [v for v in fixes.values() if isinstance(v, list)]
                fixes = extracted_list[0] if extracted_list else [fixes]

    except FileNotFoundError:
        print(f"Error: JSON file not found at {json_filepath}")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON: {e}")
        sys.exit(1)

    conn = sqlite3.connect(out_db_filepath)
    cursor = conn.cursor()

    try:
        for entry in fixes:
            # Failsafe: Ensure the entry is actually a dictionary before calling .get()
            if not isinstance(entry, dict):
                print(f"Warning: Skipping invalid entry. Expected dict, got {type(entry).__name__}: {entry}")
                continue

            entity_type = entry.get("entity_type")
            entity_id = entry.get("entity_id")
            timestamp = entry.get("timestamp")
            note = entry.get("note")
            fix_data = entry.get("fix", {})

            cursor.execute("""
                INSERT INTO manual_checks (entity_type, entity_id, timestamp, note)
                VALUES (?, ?, ?, ?)
            """, (entity_type, entity_id, timestamp, note))

            if entity_type == "compound":
                if "formula" in fix_data:
                    cursor.execute("""
                        UPDATE compounds
                        SET formula = ?, manual_flag = 1
                        WHERE compound_id = ?
                    """, (fix_data["formula"], entity_id))

            elif entity_type == "reaction":
                if "equation" in fix_data:
                    eq = fix_data["equation"]
                    
                    cursor.execute("""
                        UPDATE reactions
                        SET equation = ?, manual_flag = 1
                        WHERE reaction_id = ?
                    """, (eq, entity_id))
                    
                    substrates, products = parse_equation(eq)
                    
                    stio_dict = {}
                    for comp_id, coeff in substrates:
                        stio_dict[comp_id] = stio_dict.get(comp_id, 0.0) - coeff
                    for comp_id, coeff in products:
                        stio_dict[comp_id] = stio_dict.get(comp_id, 0.0) + coeff
                        
                    cursor.execute("""
                        DELETE FROM reaction_stoichiometry
                        WHERE reaction_id = ?
                    """, (entity_id,))
                    
                    for comp_id, net_coeff in stio_dict.items():
                        if net_coeff == 0:
                            continue # Ignore compounds that cancel out on both sides
                            
                        cursor.execute("""
                            INSERT OR IGNORE INTO compounds (compound_id, manual_flag)
                            VALUES (?, 1)
                        """, (comp_id,))
                        
                        cursor.execute("""
                            INSERT INTO reaction_stoichiometry (reaction_id, compound_id, coeff)
                            VALUES (?, ?, ?)
                        """, (entity_id, comp_id, net_coeff))
                        
                    cursor.execute("""
                        UPDATE reactions
                        SET stoichiometry = ?
                        WHERE reaction_id = ?
                    """, (json.dumps(stio_dict), entity_id))

            elif entity_type == "ec_number":
                if not fix_data:
                    continue

                if "all_reac" in fix_data and fix_data["all_reac"]:
                    cursor.execute("""
                        UPDATE ec_numbers
                        SET all_reac = ?, manual_check = 1
                        WHERE ec_number = ?
                    """, (fix_data["all_reac"], entity_id))

                # Failsafe: Safely get missing compounds and default to empty dict/lists if None
                missing_compounds = fix_data.get("missing_compounds") or {}
                substrates = missing_compounds.get("substrates") or []
                products = missing_compounds.get("products") or []

                if substrates or products:
                    for comp in substrates + products:
                        c_id = comp.get("id")
                        c_name = comp.get("name")
                        if c_id:
                            cursor.execute("""
                                INSERT OR IGNORE INTO compounds (compound_id, name, manual_flag)
                                VALUES (?, ?, 1)
                            """, (c_id, c_name))

                    cursor.execute("""
                        SELECT rowid, substrate, product 
                        FROM ec_numbers 
                        WHERE ec_number = ?
                    """, (entity_id,))
                    
                    rows = cursor.fetchall()

                    # Format missing items as readable strings, skipping entries without IDs
                    valid_subs = [f"{s.get('name', 'Unknown')} [{s.get('id')}]" for s in substrates if s.get('id')]
                    valid_prods = [f"{p.get('name', 'Unknown')} [{p.get('id')}]" for p in products if p.get('id')]
                    
                    new_subs_str = " + ".join(valid_subs)
                    new_prods_str = " + ".join(valid_prods)

                    for rowid, current_sub, current_prod in rows:
                        updated_sub = current_sub or ""  # Failsafe: handle None from DB
                        if new_subs_str:
                            updated_sub = f"{updated_sub} + {new_subs_str}" if updated_sub else new_subs_str
                            
                        updated_prod = current_prod or "" # Failsafe: handle None from DB
                        if new_prods_str:
                            updated_prod = f"{updated_prod} + {new_prods_str}" if updated_prod else new_prods_str

                        cursor.execute("""
                            UPDATE ec_numbers
                            SET substrate = ?, product = ?, manual_check = 1
                            WHERE rowid = ?
                        """, (updated_sub, updated_prod, rowid))

        # Commit all changes and close connection
        conn.commit()
        print(f"Successfully processed {len(fixes)} manual fixes.")

    except sqlite3.Error as e:
        conn.rollback() # Undo changes if something failed
        print(f"A database error occurred: {e}")
        sys.exit(1)
    finally:
        conn.close()
        
    return out_db_filepath

def main():
    parser = argparse.ArgumentParser(description="Apply manual fixes to KEGG database.")
    parser.add_argument("--kegg_db", required=True, help="Input SQLite database file")
    parser.add_argument("--manual_inst", required=True, help="JSON file with manual fixes")
    parser.add_argument("--fixed_db", required=True, help="Output SQLite database file")
    args = parser.parse_args()
    
    fixed_db = read_manual_fix(
        json_filepath=args.manual_inst,
        in_db_filepath=args.kegg_db,
        out_db_filepath=args.fixed_db
    )
    print(f"Fixed database saved to: {fixed_db}")
    
if __name__ == "__main__":
    main()