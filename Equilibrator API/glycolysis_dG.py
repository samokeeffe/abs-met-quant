import numpy as np
from equilibrator_api import ComponentContribution, Reaction, compatibility, Q_
import cobra
import re
import os
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

cc = ComponentContribution()

# Optional: changing the aqueous environment parameters
cc.p_h = Q_(7.2)
cc.p_mg = Q_(3.0)
cc.ionic_strength = Q_("0.15M")


# Create compound dictionary
compound_ids = ["h2o", "adp", "atp", "g6p", "pi", "f6p", "fdp", "nad", "3pg", "nadh", "h", "pep", "pyr", "glc__D"]
compound_dict = {cid: cc.get_compound(f"bigg.metabolite:{cid}") for cid in compound_ids}

# Define compound abundances for each cell type
abundance_cd4 = {
    "glc__D": Q_("0.77 mM"),
    "atp": Q_("8.02 mM"),
    "adp": Q_("0.413 mM"),
    "pi": Q_("3.54 mM"),
    "f6p": Q_("0.206 mM"),
    "fdp": Q_("75.3 uM"),
    "nad": Q_("0.313 mM"),
    "nah": Q_("11.8 uM"),
    "3pg": Q_("12.8 uM"),
    "pep": Q_("5.96 uM"),
    "pyr": Q_("6.15 mM"),
    "g6p": Q_("0.64 mM")
}

abundance_cd8 = {
    "glc__D": Q_("0.77 mM"),
    "atp": Q_("6.17 mM"),
    "adp": Q_("0.346 mM"),
    "pi": Q_("3.54 mM"),
    "f6p": Q_("0.151 mM"),
    "fdp": Q_("0.115 mM"),
    "nad": Q_("0.499 mM"),
    "nadh": Q_("25.2 uM"),
    "3pg": Q_("9.45 uM"),
    "pep": Q_("6.91 uM"),
    "pyr": Q_("10.1 mM"),
    "g6p": Q_("0.474 mM")
}

abundance_jurkat = {
    "glc__D": Q_("0.77 mM"),
    "atp": Q_("5 mM"),
    "adp": Q_("0.277 mM"),
    "pi": Q_("3.54 mM"),
    "f6p": Q_("0.0826 mM"),
    "fdp": Q_("113 uM"),
    "nad": Q_("0.689 mM"),
    "nadh": Q_("44.2 uM"),
    "3pg": Q_("34 uM"),
    "pep": Q_("8.49 uM"),
    "pyr": Q_("5.92 mM"),
    "g6p": Q_("0.257 mM")
}

# Map cell types to their abundances
cell_type_abundances = {
    "CD4": abundance_cd4,
    "CD8": abundance_cd8,
    "Jurkat": abundance_jurkat,
}

# Create reactions with appropriate stoichiometry

hk_reaction = Reaction({
    compound_dict["glc__D"]: -1,
    compound_dict["atp"]: -1,
    compound_dict["g6p"]: 1,
    compound_dict["adp"]: 1,
})


pgi_reaction = Reaction({
    compound_dict["g6p"]: -1,
    compound_dict["f6p"]: 1,
})

pfk_reaction = Reaction({
    compound_dict["f6p"]: -1,
    compound_dict["atp"]: -1,
    # compound_dict["h"]: 1,
    compound_dict["adp"]: 1,
    compound_dict["fdp"]: 1,
})

fba_tpi_gapdh_pgk_reaction = Reaction({
    compound_dict["fdp"]: -1,
    compound_dict["nad"]: -2,
    compound_dict["pi"]: -2,
    compound_dict["adp"]: -2,
    compound_dict["3pg"]: 2,
    compound_dict["nadh"]: 2,
    # compound_dict["h"]: 2,
    compound_dict["atp"]: 2,
})

pgam_eno_reaction = Reaction({
    compound_dict["3pg"]: -1,
    compound_dict["pep"]: 1,
    compound_dict["h2o"]: 1
})

pyk_reaction = Reaction({
    compound_dict["pep"]: -1,
    compound_dict["adp"]: -1,
    # compound_dict["h"]: -1,
    compound_dict["pyr"]: 1,
    compound_dict["atp"]: 1,
})

lumped_reaction = Reaction({
    compound_dict["glc__D"]: -1,
    compound_dict["nad"]: -2,
    compound_dict["adp"]: -2,
    compound_dict["pi"]: -2,
    compound_dict["pyr"]: 2,
    compound_dict["nadh"]: 2,
    compound_dict["atp"]: 2,
    compound_dict["h2o"]: 2,
})
# List of reactions with names
reactions = [
    ("HK", hk_reaction),
    ("PGI", pgi_reaction),
    ("PFK", pfk_reaction),
    ("FBA_TPI_GAPDH_PGK", fba_tpi_gapdh_pgk_reaction),
    ("PGAM_ENO", pgam_eno_reaction),
    ("PYK", pyk_reaction),
    ("lumped", lumped_reaction)
]

# Calculate ΔG' for each cell type and write results to a file

# Get the directory of the script being run
script_dir = os.path.dirname(os.path.abspath(__file__))

# Create the output directory next to the script
output_dir = os.path.join(script_dir, "data")
os.makedirs(output_dir, exist_ok=True)

# Set full path to the output file
output_path = os.path.join(output_dir, "deltaG_glycolysis_Tcells.txt")

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
