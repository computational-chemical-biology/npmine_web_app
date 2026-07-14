import os
from time import time
import requests
from bs4 import BeautifulSoup
import json
from rdkit import Chem
from rdkit.Chem import Draw
from flask import current_app

def cpd2prop(inchikey):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/InChIKey/%s/property/MolecularFormula,XLogP,ExactMass,Charge/JSON" % inchikey
    resp = requests.get(url)
    if resp.status_code!=200:
        d = {}
    else:
        d = json.loads(resp.text)
        d = d['PropertyTable']['Properties'][0]

    curl = 'http://classyfire.wishartlab.com/entities/%s.json' % inchikey
    resp = requests.get(curl)
    if resp.status_code!=200:
        cd = {}
    else:
        cd = json.loads(resp.text)
        if 'kingdom' in cd.keys():
            d['kingdom'] = cd['kingdom']['name']
        if 'superclass' in cd.keys():
            d['superclass'] = cd['superclass']['name']
        if 'class' in cd.keys():
            d['class'] = cd['class']['name']
        if 'subclass' in cd.keys():
            d['subclass'] = cd['subclass']['name']

    return d

def save_compound_image(compound_id, smiles):
    filename = f"{compound_id}.png"
    relative_path = os.path.join('compound_images', filename)
    
    relative_path = os.path.normpath(relative_path).replace("\\", "/")

    full_path = os.path.join(current_app.root_path, 'static', relative_path)

    os.makedirs(os.path.dirname(full_path), exist_ok=True)

    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol, size=(200, 200))
    img.save(full_path)

    return relative_path


def classify_smiles(smiles: str, logger, retries: int = 3, timeout: int = 20):

    if not smiles:
        return {
            "success": False,
            "error": "Empty SMILES",
            "data": None
        }

    last_error = None

    for attempt in range(1, retries + 1):
        try:
            response = requests.get(
                "https://npclassifier.gnps2.org/classify",
                params={"smiles": smiles},
                timeout=timeout
            )

            logger.info(
                "NPClassifier attempt=%s status=%s smiles=%s",
                attempt,
                response.status_code,
                smiles
            )

            response.raise_for_status()

            data = response.json()

            logger.info(
                "NPClassifier response=%s",
                data
            )

            return {
                "success": True,
                "error": None,
                "data": {
                    "class_results": ', '.join(
                        data.get("class_results") or []
                    ),
                    "superclass_results": ', '.join(
                        data.get("superclass_results") or []
                    ),
                    "pathway_results": ', '.join(
                        data.get("pathway_results") or []
                    ),
                    "isglycoside": data.get("isglycoside", False)
                }
            }

        except requests.Timeout as e:
            last_error = f"Timeout on attempt {attempt}"

            logger.warning(
                "NPClassifier timeout attempt=%s smiles=%s",
                attempt,
                smiles
            )

        except requests.RequestException as e:
            last_error = f"Request error: {str(e)}"

            logger.warning(
                "NPClassifier request error attempt=%s smiles=%s error=%s",
                attempt,
                smiles,
                str(e)
            )

        except ValueError as e:
            last_error = f"Invalid JSON response: {str(e)}"

            logger.warning(
                "NPClassifier invalid JSON attempt=%s smiles=%s",
                attempt,
                smiles
            )

        except Exception as e:
            last_error = f"Unexpected error: {str(e)}"

            logger.exception(
                "Unexpected NPClassifier error attempt=%s smiles=%s",
                attempt,
                smiles
            )

        if attempt < retries:
            time.sleep(attempt)

    return {
        "success": False,
        "error": last_error,
        "data": None
    }