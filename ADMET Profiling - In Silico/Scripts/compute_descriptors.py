# scripts/compute_descriptors.py

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Crippen
from rdkit.Chem import rdMolDescriptors

def compute_descriptors(smiles: str) -> dict:
    """
    Computes molecular descriptors for a given SMILES string.
    Returns a dictionary of values or None if SMILES is invalid.
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    descriptors = {
        "MW": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "Rotatable_Bonds": Lipinski.NumRotatableBonds(mol)
    }

    return descriptors


def lipinski_rule_of_five(desc: dict) -> bool:
    """
    Applies Lipinski's Rule of Five.
    Returns True if molecule passes, False otherwise.
    """

    return (
        desc["MW"] <= 500
        and desc["LogP"] <= 5
        and desc["HBD"] <= 5
        and desc["HBA"] <= 10
    )

