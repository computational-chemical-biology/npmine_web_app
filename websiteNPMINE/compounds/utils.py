import os
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