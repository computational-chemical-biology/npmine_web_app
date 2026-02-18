import requests
from websiteNPMINE.compounds.utils import save_compound_image
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
from rdkit.Chem.inchi import MolToInchi, MolToInchiKey
import pubchempy as pcp


class CompoundService:

    @staticmethod
    def create_from_smiles(smiles: str, doi_obj, user_id: int, current_app, db):
        from websiteNPMINE.models import Compounds

        smiles = (smiles or "").strip()
        if not smiles:
            return None, "SMILES vazio"

        try:
            np_response = requests.get(
                f"https://npclassifier.gnps2.org/classify",
                params={"smiles": smiles},
                timeout=20
            )

            if np_response.status_code != 200:
                return None, "NP Classifier request failed"

            np_data = np_response.json()

            class_results = ', '.join(np_data.get('class_results', []))
            superclass_results = ', '.join(np_data.get('superclass_results', []))
            pathway_results = ', '.join(np_data.get('pathway_results', []))
            isglycoside = np_data.get('isglycoside', False)

        except Exception as e:
            current_app.logger.exception("NP classifier error")
            return None, f"NP classifier error: {str(e)}"

        data, err = CompoundService.extract_from_smiles(smiles=smiles)
        
        try:
            from websiteNPMINE.compounds.routes import fetch_pubchem_data
            result = fetch_pubchem_data(inchikey=data["inchikey"])
            compound = Compounds(
                journal=None,
                compound_name=result["compound_name"],
                smiles=smiles,
                article_url=doi_obj.doi,
                inchi_key=data["inchikey"],
                exact_molecular_weight=data["exact_molecular_weight"],
                class_results=class_results,
                superclass_results=superclass_results,
                pathway_results=pathway_results,
                isglycoside=isglycoside,
                pubchem_id=result["pubchem_id"],
                inchi=data["inchi"],
                source='NPMine',
                user_id=user_id
            )

            db.session.add(compound)
            db.session.commit()

            compound.dois.append(doi_obj)
            db.session.commit()

            img_path = save_compound_image(compound.id, smiles)
            if img_path:
                compound.compound_image = img_path
                db.session.commit()

            return compound, None

        except Exception as e:
            db.session.rollback()
            current_app.logger.exception("Error creating compound from SMILES")
            return None, f"DB error: {str(e)}"
    
    @staticmethod
    def extract_from_smiles(smiles: str):

        if not smiles:
            return None, "SMILES vazio"

        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None, "RDKit não conseguiu parsear SMILES"

            smiles_canonical = Chem.MolToSmiles(mol, canonical=True)

            inchi = MolToInchi(mol)
            inchikey = MolToInchiKey(mol)

            exact_mw = Descriptors.ExactMolWt(mol)
            formula = rdMolDescriptors.CalcMolFormula(mol)
            logp = Descriptors.MolLogP(mol)
            num_atoms = mol.GetNumAtoms()
            ring_count = rdMolDescriptors.CalcNumRings(mol)

            fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

            result = {
                "smiles": smiles_canonical,
                "inchi": inchi,
                "inchikey": inchikey,
                "exact_molecular_weight": exact_mw,
                "molecular_formula": formula,
                "logp": logp,
                "num_atoms": num_atoms,
                "ring_count": ring_count,
                "mol": mol,
                "fingerprint": fingerprint
            }

            return result, None

        except Exception as e:
            return None, f"RDKit error: {str(e)}"

    @staticmethod
    def smiles_from_inchikey(inchikey: str):

        if not inchikey:
            return None, "InChIKey vazio"

        try:
            compounds = pcp.get_compounds(inchikey, namespace="inchikey")

            if not compounds:
                return None, "Nenhum composto encontrado no PubChem"

            comp = compounds[0]

            smiles = comp.canonical_smiles
            inchi = comp.inchi
            cid = comp.cid

            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None, "SMILES retornado inválido no RDKit"

            return {
                "smiles": smiles,
                "inchi": inchi,
                "pubchem_id": cid
            }, None

        except Exception as e:
            return None, f"Erro resolvendo InChIKey: {str(e)}"