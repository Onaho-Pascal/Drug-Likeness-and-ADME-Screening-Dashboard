# scripts/admet_prediction.py

import requests

def get_admet_predictions(smiles: str) -> dict:
    """
    Sends a SMILES string to ADMETlab 2.0 API and retrieves:
    - Human Intestinal Absorption (HIA)
    - Blood-Brain Barrier Penetration (BBB)
    - Hepatotoxicity
    - hERG inhibition

    Returns a dictionary, or None if request fails.
    """

    url = "https://admet.scbdd.com/predict/"
    payload = {"smiles": smiles}

    try:
        response = requests.post(url, json=payload, timeout=10)

        if response.status_code != 200:
            return None

        data = response.json()

        # Extract key ADMET values
        return {
            "HIA": data.get("hia", None),                     # Absorption
            "BBB": data.get("bbb", None),                     # Brain permeability
            "Hepatotoxicity": data.get("hepatotoxicity", None),
            "hERG_Block": data.get("herg", None)
        }

    except Exception as e:
        print(f"ADMET API Error: {e}")
        return None

