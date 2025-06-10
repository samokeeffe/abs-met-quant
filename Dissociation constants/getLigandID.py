from zeep import Client
import hashlib
import csv

# Define WSDL URL
wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"

# Placeholder email and password 
username = "j.doe@example.edu"
password = hashlib.sha256("password".encode("utf-8")).hexdigest()

# Initialize the SOAP client
client = Client(wsdl)


# List of metabolites to search, select metabolites of interest
metabolites = [
    "1,3-bisphosphoglycerate",
    "2,3-bisphosphoglycerate",
    "2,3-dihydroxybenzoic acid",
    "2-dehydro-D-gluconate",
    "2-phosphoglycerate",
    "3-phosphoglycerate",
    "3-phospho-serine",
    "4-hydroxybenzoate",
    "6-phospho-D-gluconate",
    "acetoacetyl-CoA",
    "acetyl-CoA",
    "acetylphosphate",
    "aconitate",
    "ac-serine",
    "adenine",
    "adenosine",
    "adenosine-phosphosulfate",
    "ADP",
    "ADP-glucose",
    "a-ketoglutarate",
    "alanine",
    "AMP",
    "anthranilate",
    "arginine",
    "asparagine",
    "aspartate",
    "ATP",
    "carbamoyl-aspartate",
    "citrate",
    "citrulline",
    "CMP",
    "coenzyme-A",
    "CTP",
    "cyclic-AMP",
    "cysteine",
    "cytidine",
    "cytosine",
    "dAMP",
    "dATP",
    "dCDP",
    "dCMP",
    "dCTP",
    "deoxyadenosine",
    "deoxyguanosine",
    "deoxyribose-5-phosphate",
    "dGMP",
    "dihydroorotate",
    "dihydroxyacetonephosphate",
    "dTDP",
    "dTMP",
    "dTTP",
    "erythrose-4-phosphate",
    "FAD",
    "flavin mononucleotide",
    "fructose-1,6-bisphosphate",
    "fructose-6-phosphate",
    "fumarate",
    "GDP",
    "gluconate",
    "gluconolactone",
    "glucosamine-6-phosphate",
    "glucose-6-phosphate",
    "glutamate",
    "glutamine",
    "glutathione",
    "glutathione disulfide",
    "glyceraldehyde 3-phosphate",
    "glycerate",
    "glycine",
    "GMP",
    "GTP",
    "guanine",
    "guanosine",
    "histidine",
    "histidinol",
    "homocysteine",
    "hydroxyisocaproic acid",
    "IDP",
    "IMP",
    "inosine",
    "isocitrate",
    "isoleucine",
    "ITP",
    "leucine",
    "lysine",
    "malate",
    "malonyl-CoA",
    "methionine",
    "methylmalonic acid",
    "myo-inositol",
    "N-acetyl-glucosamine-1/6-phosphate",
    "N-acetyl-glutamine",
    "N-Acetyl-L-alanine",
    "N-Acetyl-L-aspartic acid",
    "N-acetyl-ornithine",
    "NAD+",
    "NADH",
    "NADP+",
    "NADPH",
    "ornithine",
    "orotate",
    "oxaloacetate",
    "phenylalanine",
    "phenylpyruvate",
    "phosphoenolpyruvate",
    "proline",
    "propionyl-CoA",
    "PRPP",
    "pyruvate",
    "quinolinate",
    "riboflavin",
    "ribose-5-phosphate",
    "ribulose-5-phosphate",
    "S-adenosyl-L-homocysteine",
    "S-adenosyl-L-methionine",
    "sedoheptulose-7-phosphate",
    "serine",
    "shikimate",
    "sn-glycerol 3-phosphate",
    "succinate",
    "succinyl-CoA",
    "taurine",
    "threonine",
    "thymidine",
    "trehalose",
    "tryptophan",
    "tyrosine",
    "UDP",
    "UDP-glucose",
    "UDP-glucuronate",
    "UDP-N-acetyl-glucosamine",
    "UMP",
    "uridine",
    "UTP",
    "valine",
    "xylulose-5-phosphate"
    "D-glucosamine 6-phosphate"
]

# Dictionary to store results
ligand_ids = {}

# Loop through each metabolite and fetch ligandStructureId
for metabolite in metabolites:
    try:
        parameters = (username, password, metabolite)
        ligand_id = client.service.getLigandStructureIdByCompoundName(*parameters)
        
        if ligand_id:
            ligand_ids[metabolite] = ligand_id
            print(f"{metabolite}: {ligand_id}")
        else:
            print(f"{metabolite}: No ID found")
    except Exception as e:
        print(f"Error fetching {metabolite}: {e}")

