
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import rdMolDescriptors

data = {
    "Drug": [
        "Chloroquine", "Hydroxychloroquine", 
        "Amodiaquine", "Artemisinin", "Lumefantrine"
    ],
    "SMILES": [
        "CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12",
        "CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12O",
        "CCNc1ccc(OCCN(CC)CC)nc2c(Cl)cccc12",
        "CC1CCC2C(C1)C(OO2)C(=O)O",
        "CC(C)(C)OC1=C(C=C(C=C1)OC2CCCC2)OCCN3CCN(CC3)C"
    ]
}

df = pd.DataFrame(data)
df

def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    return {
        "MW": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "Rotatable_bonds": Lipinski.NumRotatableBonds(mol)
    }

def lipinski_pass(desc):
    return (
        desc["MW"] <= 500 and
        desc["LogP"] <= 5 and
        desc["HBD"] <= 5 and
        desc["HBA"] <= 10
    )
results = []
for i, row in df.iterrows():
    desc = compute_descriptors(row["SMILES"])
    if desc:
        desc["Lipinski"] = lipinski_pass(desc)
        desc["Drug"] = row["Drug"]
        results.append(desc)

results_df = pd.DataFrame(results)
results_df

results_df.to_csv("druglikeness_results.csv", index=False)

import requests

def admet_prediction(smiles):
    url = "https://admet.scbdd.com/predict/"
    payload = {"smiles": smiles}
    
    try:
        response = requests.post(url, json=payload)
        data = response.json()
        return {
            "Human_Intestinal_Absorption": data["hia"],
            "BBB_Penetration": data["bbb"],
            "Hepatotoxicity": data["hepatotoxicity"],
            "hERG_Block": data["herg"]
        }
    except:
        return None

# Apply to dataset
admet_list = []
for i, row in df.iterrows():
    pred = admet_prediction(row["SMILES"])
    if pred:
        pred["Drug"] = row["Drug"]
        admet_list.append(pred)

admet_df = pd.DataFrame(admet_list)
admet_df

final_df = results_df.merge(admet_df, on="Drug", how="left")
final_df

final_df.to_csv("drug_chemoinformatics_summary.csv", index=False)
