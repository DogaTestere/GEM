import cobra
from cobra.io import read_sbml_model, write_sbml_model
import os
import argparse

def analyze_and_balance_model(xml_path):
    """
    Reads a boundary-free cobraPy model, tallies stoichiometric coefficients, 
    adds appropriate bounds/exchanges, and returns the balanced model and missing formulas.
    """
    
    if not os.path.exists(xml_path):
        raise FileNotFoundError(f"Critical Error: XML file '{xml_path}' not found. Crashing.")

    model = read_sbml_model(xml_path)
    
    missing_formulas = set()
    
    for met in model.metabolites:
        if not met.formula:
            missing_formulas.add(met.id)
            
        # Tally the stoichiometric coefficients across all reactions
        # e.g., (+2) + (-2) + (+1) + (+1) = +2
        tally = sum(rxn.get_coefficient(met) for rxn in met.reactions)
        
        if tally > 0:
            # Overflow -> Add a Sink to consume the excess
            model.add_boundary(met, type="sink", lb=-tally, ub=0)
            
        elif tally < 0:
            # Deficiency -> Add a Sink to supply the missing amount
            model.add_boundary(met, type="sink", lb=0, ub=-tally)

    return model, missing_formulas

def main(xml_path, output_xml_path, report_path):
    try:
        model, missing_formulas = analyze_and_balance_model(xml_path)
    except FileNotFoundError as e:
        print(e)
        return

    with open(report_path, "w", encoding="utf-8") as f:
        if missing_formulas:
            f.write("Metabolites missing formulas:\n")
            for met_id in sorted(missing_formulas):
                f.write(f"{met_id}\n")
        else:
            f.write("All metabolites have chemical formulas.\n")
            
    print(f"\nMissing formulas report saved to: {report_path}")
    
    write_sbml_model(model, output_xml_path)
    print(f"Balanced model saved to: {output_xml_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--input_model", required=True)
    parser.add_argument("--output_model", required=True)
    parser.add_argument("--report_out")
    
    args = parser.parse_args()
    
    main(
        xml_path=args.input_model,
        output_xml_path=args.output_model,
        report_path=args.report_out
    )