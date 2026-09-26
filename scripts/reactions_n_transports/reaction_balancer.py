#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sqlite3
import os
import sys
import argparse
import shutil
import warnings
import logging
import json 

# Uyarıları tamamen kapat
warnings.filterwarnings("ignore")

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(levelname)s: %(message)s',
        stream=sys.stdout
    )

def run_thermo_update(input_db, output_db, missing_json_path):
    if not os.path.exists(input_db):
        logging.error(f"Giriş veritabanı bulunamadı: {input_db}")
        sys.exit(1)

    # Çıktı veritabanını oluştur
    shutil.copy2(input_db, output_db)
    conn = sqlite3.connect(output_db)
    cursor = conn.cursor()
    
    logging.info("eQuilibrator API başlatılıyor...")
    try:
        from equilibrator_api import ComponentContribution
        cc = ComponentContribution()
    except Exception as e:
        logging.error(f"API başlatılamadı: {e}")
        sys.exit(1)

    # Gerekli sütunları ekle
    try:
        cursor.execute("ALTER TABLE reactions ADD COLUMN lower_bound FLOAT")
        cursor.execute("ALTER TABLE reactions ADD COLUMN upper_bound FLOAT")
    except sqlite3.OperationalError:
        pass

    cursor.execute("SELECT reaction_id FROM reactions")
    reactions = cursor.fetchall()
    
    processed_count = 0
    skipped_count = 0
    
    error_report = {
        "summary": {},
        "details": []
    }

    logging.info(f"Toplam {len(reactions)} reaksiyon işleniyor...")

    for (reaction_id,) in reactions:
        cursor.execute(
            "SELECT compound_id, coeff FROM reaction_stoichiometry WHERE reaction_id=?",
            (reaction_id,)
        )
        stoich_rows = cursor.fetchall()

        left, right = [], []
        can_process = True
        problematic_compounds = []
        
        for compound_id, coeff in stoich_rows:
            kegg_id = f"kegg:{compound_id}"
            try:
                comp = cc.get_compound(kegg_id)
                if comp is None:
                    raise ValueError(f"Cache miss")
                
                val = float(coeff)
                if val < 0:
                    left.append(f"{abs(val)} {kegg_id}")
                elif val > 0:
                    right.append(f"{abs(val)} {kegg_id}")
            except Exception:
                problematic_compounds.append(compound_id)
                can_process = False

        lower, upper = -1000.0, 1000.0
        success = False

        if can_process and left and right:
            try:
                reaction_string = " + ".join(left) + " = " + " + ".join(right)
                rxn = cc.parse_reaction_formula(reaction_string)
                dG_prime = cc.standard_dg_prime(rxn)
                
                mag = dG_prime.to("kJ/mol").magnitude
                dG_value = float(mag.nominal_value) if hasattr(mag, "nominal_value") else float(mag)

                if dG_value < -3.0:
                    lower, upper = 0.0, 1000.0
                elif dG_value > 3.0:
                    lower, upper = -1000.0, 0.0
                else:
                    lower, upper = -1000.0, 1000.0
                
                success = True
                processed_count += 1
            except Exception as e:
                error_report["details"].append({
                    "reaction_id": reaction_id,
                    "reason": "Calculation error",
                    "error_message": str(e)
                })
        else:
            if problematic_compounds:
                error_report["details"].append({
                    "reaction_id": reaction_id,
                    "reason": "Missing compounds",
                    "compounds": problematic_compounds
                })
            else:
                error_report["details"].append({
                    "reaction_id": reaction_id,
                    "reason": "Empty stoichiometry/Imbalanced"
                })

        cursor.execute(
            "UPDATE reactions SET lower_bound=?, upper_bound=? WHERE reaction_id=?",
            (lower, upper, reaction_id)
        )
        
        if not success:
            skipped_count += 1

    conn.commit()
    conn.close()

    error_report["summary"] = {
        "total_reactions": len(reactions),
        "successfully_calculated": processed_count,
        "default_bounds_assigned": skipped_count
    }

    # JSON raporunu kaydet
    with open(missing_json_path, "w", encoding="utf-8") as jf:
        json.dump(error_report, jf, indent=4, ensure_ascii=False)

    logging.info(f"JSON raporu '{missing_json_path}' kaydedildi.")
    logging.info(f"İşlem Tamamlandı: {processed_count} başarılı, {skipped_count} varsayılan.")
    
    sys.exit(0)

def main():
    setup_logging()
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_db", required=True)
    parser.add_argument("--balanced_db", required=True)
    parser.add_argument("--missing_json", required=True) # Argüman burada tanımlı
    args = parser.parse_args()
    
    # Argümanı fonksiyona paslıyoruz:
    run_thermo_update(args.input_db, args.balanced_db, args.missing_json)

if __name__ == "__main__":
    main()