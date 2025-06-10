import numpy as np
from equilibrator_api import ComponentContribution, Reaction, compatibility, Q_
import cobra
import re
import os

cc = ComponentContribution()

os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Optional: changing the aqueous environment parameters
cc.p_h = Q_(7.2)
cc.p_mg = Q_(3.0)
cc.ionic_strength = Q_("0.15M")


# Create compound dictionary
compound_ids = ["nadph", "nadp", "gthox", "gthrd", "h"]
compound_dict = {cid: cc.get_compound(f"bigg.metabolite:{cid}") for cid in compound_ids}

# Define compound abundances for each cell type
abundance_cd4 = {
    "nadph": Q_("4.41 uM"),
    "nadp" : Q_("0.2 uM"),
    "gthox": Q_("0.63 mM"),
    "gthrd" : Q_("5.63 mM")
}

abundance_cd8 = {
    "nadph": Q_("5.18 uM"),
    "nadp" : Q_("0.33 uM"),
    "gthox": Q_("0.83 mM"),
    "gthrd" : Q_("7.08 mM")
}

abundance_jurkat = {
    "nadph": Q_("3.59 uM"),
    "nadp" : Q_("0.16 uM"),
    "gthox": Q_("0.91 mM"),
    "gthrd" : Q_("6.86 mM")
}

# Map cell types to their abundances
cell_type_abundances = {
    "CD4": abundance_cd4,
    "CD8": abundance_cd8,
    "Jurkat": abundance_jurkat,
}

# Create reactions with appropriate stoichiometry

glutathione_reductase = Reaction({
    compound_dict["nadph"]: -1,
    compound_dict["gthox"]: -1,
    compound_dict["nadp"]: 1,
    compound_dict["gthrd"]: 2,
    compound_dict["h"]: -1
})

# List of reactions with names
reactions = [
    ("Glutathione Reductase", glutathione_reductase)
]

# Calculate ΔG' for each cell type and write results to a file

# Get the directory of the script being run
script_dir = os.path.dirname(os.path.abspath(__file__))

# Create the output directory next to the script
output_dir = os.path.join(script_dir, "data")
os.makedirs(output_dir, exist_ok=True)

# Set full path to the output file
output_path = os.path.join(output_dir, "deltaG_glutathione_Tcells.txt")
with open(output_path, 'w') as outF:
    for cell_type, abundances in cell_type_abundances.items():
        for name, reaction in reactions:
            # Apply abundances for the current cell type
            for compound, abundance in abundances.items():
                if compound in compound_dict:
                    reaction.set_abundance(compound_dict[compound], abundance)
            # Calculate ΔG'
            dG_prime = cc.dg_prime(reaction)
            deltaG_value = dG_prime.value.m_as("kJ/mol")
            deltaG_error = dG_prime.error.m_as("kJ/mol")
            output_s = (f"Cell Type={cell_type} ID={name}  "
                        f"ΔG'={deltaG_value:.2f} kJ/mol  "
                        f"ΔG'_error={deltaG_error:.2f} kJ/mol\n")
            print(output_s)
            outF.write(output_s)
