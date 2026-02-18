from flask import Blueprint, render_template, flash, redirect, url_for, request, current_app, g,jsonify, Response, abort
from flask_login import login_required, current_user
from websiteNPMINE.compounds.compound_service import CompoundService
from websiteNPMINE.compounds.forms import CompoundForm, SearchForm, CompoundEditForm
from websiteNPMINE.models import Compounds,DOI,Taxa, CompoundHistory
from websiteNPMINE import db
import os
import json
import requests
from websiteNPMINE.compounds.utils import *
from rdkit import Chem
from rdkit.Chem import Draw, DataStructs, AllChem
import pubchempy as pcp
from PIL import Image
from datetime import datetime, timezone
from flask_babel import Babel, get_locale
from rdkit.Chem.Fingerprints import FingerprintMols
from werkzeug.utils import secure_filename
from urllib.parse import quote
from websiteNPMINE import csrf
from io import StringIO
import csv

compounds = Blueprint('compounds', __name__)

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

import requests

@compounds.route('/new_compound', methods=['GET', 'POST'])
@login_required
def registerCompound():
    logged_in = current_user.is_authenticated
    form = CompoundForm()

    if form.validate_on_submit():
        doi = form.doi.data

        if not doi:
            flash('DOI Link is required', 'error')
            return redirect(request.url)

        existing_doi = DOI.query.filter_by(doi=doi).first()
        if not existing_doi:
            new_doi = DOI(doi=doi)
            db.session.add(new_doi)
            db.session.commit()
            existing_doi = new_doi
        else:   
            flash('DOI already in database!', 'info')

        compound_blocks = request.form.getlist('inchikey')
        smiles_blocks = request.form.getlist('smiles')
        for i, inchikey in enumerate(compound_blocks):
            smiles_input = smiles_blocks[i].strip() if i < len(smiles_blocks) else None

            genus = request.form.getlist('genus')[i]
            species = request.form.getlist('species')[i]
            origin_type = request.form.getlist('origin_type')[i]

            inchikey = (inchikey or "").strip()

            if not inchikey and not smiles_input:
                flash(f'Provide InChIKey or SMILES for Compound {i + 1}', 'error')
                continue

            if inchikey:
                existing_compound = Compounds.query.filter_by(inchi_key=inchikey).first()

            if not inchikey and smiles_input:
                compound, err = CompoundService.create_from_smiles(
                    smiles=smiles_input,
                    doi_obj=existing_doi,
                    user_id=current_user.id,
                    current_app=current_app,
                    db=db
                )

                if err:
                    flash(f'Compound {i+1}: {err}', 'error')
                    continue

                flash(f'Compound {i+1} created from SMILES', 'success')
                continue
            
            if existing_compound:
                # Check if the current user already has this compound
                user_compound = Compounds.query.filter_by(inchi_key=inchikey, user_id=current_user.id).first()
                if user_compound:
                    flash(f'You already have this compound with InChI Key {inchikey} in your records.', 'info')
                else:
                    # Create a new compound entry for the current user
                    new_compound = Compounds(
                        journal=existing_compound.journal,
                        compound_name=existing_compound.compound_name,
                        smiles=existing_compound.smiles,
                        article_url=existing_compound.article_url,
                        inchi_key=existing_compound.inchi_key,
                        exact_molecular_weight=existing_compound.exact_molecular_weight,
                        class_results=existing_compound.class_results,
                        superclass_results=existing_compound.superclass_results,
                        pathway_results=existing_compound.pathway_results,
                        isglycoside=existing_compound.isglycoside,
                        pubchem_id=existing_compound.pubchem_id,
                        inchi=existing_compound.inchi,
                        source=existing_compound.source,
                        user_id=current_user.id  # Change the user_id to the current user
                    )
                    db.session.add(new_compound)
                    db.session.commit()

                    # Associate the new compound with the DOI
                    new_compound.dois.append(existing_doi)
                    db.session.commit()

                    # Save the compound image for the duplicated compound
                    img_path = save_compound_image(new_compound.id, existing_compound.smiles)
                    if img_path:
                        new_compound.compound_image = img_path
                        db.session.commit()

                    flash(f'Compound with InChI Key {inchikey} duplicated for your account.', 'info')
            else:
                pubchem_data = fetch_pubchem_data(inchikey)

                if not pubchem_data:
                    flash(f'Failed to fetch data from PubChem for Compound', 'error')
                    return render_template('new_compound.html', form=form, logged_in=logged_in)

                smiles = pubchem_data.get('smiles')
                if not smiles:
                    flash(f'SMILES not found for Compound {i + 1}', 'error')
                    return render_template('new_compound.html', form=form, logged_in=logged_in)

                # Fetch NP Classifier data
                try:
                    np_response = requests.get(f"https://npclassifier.gnps2.org/classify?smiles={smiles}")
                    if np_response.status_code == 200:
                        np_data = np_response.json()
                        class_results = ', '.join(np_data.get('class_results', []))
                        superclass_results = ', '.join(np_data.get('superclass_results', []))
                        pathway_results = ', '.join(np_data.get('pathway_results', []))
                        isglycoside = np_data.get('isglycoside', False)
                    else:
                        flash(f'Failed to fetch NP Classifier data for Compound {i + 1}', 'error')
                        continue
                except Exception as e:
                    flash(f'Error fetching NP Classifier data: {e}', 'error')
                    continue

                compound = Compounds(
                    journal=None,
                    compound_name=pubchem_data['compound_name'],
                    smiles=smiles,
                    article_url=doi,
                    inchi_key=inchikey,
                    exact_molecular_weight=pubchem_data['exact_molecular_weight'],
                    class_results=class_results,
                    superclass_results=superclass_results,
                    pathway_results=pathway_results,
                    isglycoside=isglycoside,
                    pubchem_id=pubchem_data['pubchem_id'],
                    inchi=pubchem_data['inchi'],
                    source='NPMine',
                    user_id=current_user.id
                )
                db.session.add(compound)
                db.session.commit()

                compound.dois.append(existing_doi)
                db.session.commit()

                img_path = save_compound_image(compound.id, smiles)
                if img_path:
                    compound.compound_image = img_path
                    db.session.commit()

            if genus and species:
                existing_taxon = Taxa.query.filter_by(verbatim=f"{genus} {species}").first()
                if not existing_taxon:
                    taxon = Taxa(
                        article_url=doi,
                        verbatim=f"{genus} {species}",
                        user_id=current_user.id
                    )
                    db.session.add(taxon)
                    db.session.commit()

                    taxon.dois.append(existing_doi)
                else:
                    if existing_taxon not in existing_doi.taxa:
                        existing_doi.taxa.append(existing_taxon)
                
                db.session.commit()

        flash('Compounds added successfully!', 'success')
        return redirect(url_for('compounds.registerCompound'))

    for error in form.errors.values():
        flash(error[0], 'error')

    return render_template('new_compound.html', form=form, logged_in=logged_in)


def fetch_pubchem_data(inchikey):
    try:
        required_properties = ['canonical_smiles', 'inchi', 'molecular_weight', 'synonyms']
        compound = pcp.get_compounds(inchikey, 'inchikey', properties=required_properties)
        
        if not compound:
            print(f"No compound found for InChIKey: {inchikey}", flush=True)
            return None
        
        compoundData = compound[0]

        record = compoundData.to_dict(properties=[
            'canonical_smiles',
            'inchi',
            'molecular_weight',
            'iupac_name'
        ])
        
        smiles = record.get('canonical_smiles')
        inchi = record.get('inchi')
        mw = record.get('molecular_weight')
        name = record.get('iupac_name')
       
        smiles = compoundData.canonical_smiles
        name = compoundData.synonyms[0] if compoundData.synonyms else None

        if not smiles or not name:
            print(f"DEBUG: Atributos padrão falharam para {inchikey}. Iniciando extração manual...", flush=True)
            
            raw_props = vars(compoundData).get('_record', {}).get('props', [])
            
            for p in raw_props:
                label = p.get('urn', {}).get('label')
                value = p.get('value', {}).get('sval') or p.get('value', {}).get('fval')
                
                if not smiles and label == 'SMILES':
                    smiles = value
                
                if not name and label == 'IUPAC Name':
                    name = value

        return {
            'compound_name': name,
            'smiles': smiles,
            'inchi': compoundData.inchi,    
            'exact_molecular_weight': compoundData.molecular_weight,
            'pubchem_id': compoundData.cid
        }

    except Exception as e:
        print(f"Erro na extração: {str(e)}", flush=True)
        return {}

@compounds.route('/search_menu')
def search_menu():
    logged_in = current_user.is_authenticated
    return render_template('search_menu.html', logged_in=logged_in)

@compounds.route('/search')
@csrf.exempt
def search():
    logged_in = current_user.is_authenticated
    q = request.args.get("q")

    if q:
        results = Compounds.query \
            .outerjoin(Compounds.dois) \
            .filter(
                (Compounds.compound_name.ilike(f"%{q}%") |
                Compounds.smiles.ilike(f"%{q}%") |
                Compounds.inchi_key.ilike(f"%{q}%") |
                DOI.doi.ilike(f"%{q}%")) &
                (Compounds.status == 'public')  
            ) \
            .order_by(Compounds.compound_name.asc()) \
            .limit(100) \
            .all()
    else:
        results = []

    if request.headers.get('HX-Request') == 'true':
        return render_template('search_results.html', results=results)
    
    return render_template('search.html', results=results, logged_in=logged_in)

@compounds.route('/search_doi')
@csrf.exempt
def search_doi():
    logged_in = current_user.is_authenticated
    q = request.args.get("q")

    if q:
        results = (
            DOI.query
            .join(DOI.compounds)  
            .filter(Compounds.status == 'public')  
            .filter(DOI.doi.ilike(f"%{q}%")) 
            .distinct()
            .options(
                db.joinedload(DOI.compounds.and_(Compounds.status == 'public')),  
                db.joinedload(DOI.taxa)  
            )
            .order_by(DOI.doi.asc())  
            .limit(100)
            .all()
        )
    else:
        results = []


    if request.headers.get('HX-Request') == 'true':
        return render_template('search_results_doi.html', results=results)

    return render_template('search_doi.html', results=results, logged_in=logged_in)

@compounds.route('/search_taxon')
@csrf.exempt
def search_taxon():
    logged_in = current_user.is_authenticated
    q = request.args.get("q")

    if q:
      
        results = (
            Taxa.query
            .join(Taxa.dois) 
            .join(DOI.compounds)  
            .filter(Compounds.status == 'public') 
            .filter(Taxa.verbatim.ilike(f"%{q}%"))  
            .distinct()
            .options(
                db.joinedload(Taxa.dois).joinedload(
                    DOI.compounds.and_(Compounds.status == 'public')  
                ),
                db.joinedload(Taxa.dois)  
            )
            .order_by(Taxa.verbatim.asc())  
            .limit(100)
            .all()
        )
    else:
        results = []

    
    if request.headers.get('HX-Request') == 'true':
        return render_template('search_results_taxon.html', results=results)

    return render_template('search_taxon.html', results=results, logged_in=logged_in)

def pairSimilarity(inchipair, metric=DataStructs.TanimotoSimilarity, fptype='circular'):
    """ Calculates fingerprint similarity between two InChIs """
    ms = [Chem.MolFromInchi(x) for x in inchipair]
    if fptype == 'circular':
        fps = [AllChem.GetMorganFingerprint(x, 2) for x in ms]
    else:
        fps = [FingerprintMols.FingerprintMol(x) for x in ms]
    return metric(fps[0], fps[1])


def process_similarity_search(inchi, cpds):
    """Process compound similarity search given an InChI and compound list."""
    search_res = []
    for c in cpds:
        try:
            tsc = pairSimilarity([c.inchi, inchi])
        except Exception as e:
            print(f"Error calculating similarity for compound {c.compound_name}: {e}")
            continue
        
        link = f'<a href="/compound/{c.id}">{c.compound_name}</a>'  
        search_res.append([link, c.inchi_key, tsc])
    
  
    search_res.sort(key=lambda x: x[2], reverse=True)
    
    return search_res

@compounds.route('/search_structure', methods=['GET', 'POST'])
@csrf.exempt
def search_structure():
    logged_in = current_user.is_authenticated

    if request.method == 'GET':
        return render_template('search_structure.html')

  
    cpds = db.session.query(Compounds).filter(Compounds.status == 'public').all()
    search_res = []
    query = None

  
    if 'textSmiles' in request.form and request.form['textSmiles']:
        smiles = request.form['textSmiles']
        try:
            inchi = Chem.MolToInchi(Chem.MolFromSmiles(smiles))
            query = f"https://cactus.nci.nih.gov/chemical/structure/{quote(inchi)}/image"
            print(f"Generated Image URL: {query}") 
            search_res = process_similarity_search(inchi, cpds)
        except Exception as e:
            print(f"Error processing SMILES input: {e}")
            return "Invalid SMILES input", 400


    elif 'drawSmiles' in request.form and request.form['drawSmiles']:
        smiles = request.form['drawSmiles']
        try:
            inchi = Chem.MolToInchi(Chem.MolFromSmiles(smiles))
            query = f"https://cactus.nci.nih.gov/chemical/structure/{quote(inchi)}/image"
            print(f"Generated Image URL: {query}") 
            search_res = process_similarity_search(inchi, cpds)
        except Exception as e:
            print(f"Error processing drawn SMILES: {e}")
            return "Invalid SMILES input", 400

    elif 'data' in request.files:
        data = request.files['data']
        if data.filename == '':
            return "No selected file", 400

        
        try:
            file_content = data.read().decode('utf-8')
            smiles_list = file_content.splitlines()
        except Exception as e:
            print(f"Error reading file: {e}")
            return "Error reading file", 400

        
        for smiles in smiles_list:
            if smiles.strip():
                try:
                    inchi = Chem.MolToInchi(Chem.MolFromSmiles(smiles))
                    query = f"https://cactus.nci.nih.gov/chemical/structure/{quote(inchi)}/image"
                    print(f"Generated Image URL: {query}")  # Log the generated URL for debugging
                    search_res += process_similarity_search(inchi, cpds)
                except Exception as e:
                    print(f"Error processing SMILES '{smiles}': {e}")
                    continue

 
    else:
        return "Invalid request", 400

    return render_template('search_results_structure.html', dfinal=search_res, query=query, logged_in=logged_in)


def update_compound_relationships(compound, doi_data_list):
    for existing_doi in compound.dois:
        matching_doi_data = next((doi_data for doi_data in doi_data_list if doi_data['doi'] == existing_doi.doi), None)
        if matching_doi_data:
            existing_doi.some_attribute = matching_doi_data.get('some_attribute', existing_doi.some_attribute)
            doi_data_list.remove(matching_doi_data)
    
    for new_doi_data in doi_data_list:
        new_doi = DOI(**new_doi_data)
        compound.dois.append(new_doi)
    
    db.session.commit()

@compounds.route('/compound/<int:id>/edit', methods=['GET', 'POST'])
@login_required
def edit_compound(id):
    logged_in = current_user.is_authenticated
    compound = Compounds.query.get_or_404(id)

    is_admin = current_user.role_id == 1
    is_editor = current_user.role_id == 2
    is_owner = compound.user_id == current_user.id

    if not (is_admin or is_editor or is_owner):
        flash('You do not have permission to edit this compound.', 'danger')
        return redirect(url_for('main.index'))

    form = CompoundEditForm(obj=compound)

    history_records = CompoundHistory.query.filter_by(compound_id=id).order_by(CompoundHistory.created_at.desc()).all()

    form.dois.entries.clear()
    for doi in compound.dois:
        form.dois.append_entry({'doi': doi.doi})

    related_taxa = [taxa for doi in compound.dois for taxa in doi.taxa]

    if form.validate_on_submit():
        print("Form data:", form.data)

        old_smiles = compound.smiles

        history_entry = CompoundHistory(
            compound_id=compound.id,
            journal=compound.journal,
            smiles=compound.smiles,
            article_url=compound.article_url,
            inchi_key=compound.inchi_key,
            exact_molecular_weight=compound.exact_molecular_weight,
            class_results=compound.class_results,
            superclass_results=compound.superclass_results,
            pathway_results=compound.pathway_results,
            pubchem_id=compound.pubchem_id,
            inchi=compound.inchi,
        )

        db.session.add(history_entry)

        compound.journal = form.journal.data or compound.journal
        compound.compound_name = form.compound_name.data or compound.compound_name
        compound.smiles = form.smiles.data or compound.smiles
        compound.exact_molecular_weight = form.exact_molecular_weight.data or compound.exact_molecular_weight
        compound.class_results = form.class_results.data or compound.class_results
        compound.superclass_results = form.superclass_results.data or compound.superclass_results
        compound.pathway_results = form.pathway_results.data or compound.pathway_results
        compound.isglycoside = (form.isglycoside.data == 'True') 
        compound.pubchem_id = form.pubchem_id.data or compound.pubchem_id
        compound.inchi = form.inchi.data or compound.inchi
        compound.article_url = form.article_url.data or compound.article_url
        compound.ispublic = (form.ispublic.data == 'True')
        compound.inchi_key = form.inchi_key.data or compound.inchi_key

        print("DOIs before update:", [d.doi for d in compound.dois])
        
        form_doi_data_list = [entry.data['doi'] for entry in form.dois.entries]

        for existing_doi in compound.dois[:]:
            if existing_doi.doi in form_doi_data_list:
                form_doi_entry = next(entry for entry in form.dois.entries if entry.data['doi'] == existing_doi.doi)
            else:
                compound.dois.remove(existing_doi)
                db.session.delete(existing_doi)

        for form_doi_entry in form.dois.entries:
            form_doi_data = form_doi_entry.data
            if not any(doi.doi == form_doi_data['doi'] for doi in compound.dois):
                new_doi = DOI(doi=form_doi_data['doi'])
                compound.dois.append(new_doi)
                db.session.add(new_doi)

        if compound.smiles != old_smiles:
            print("SMILES changed, updating image...")
            compound_image_path = save_compound_image(compound.id, compound.smiles)
            compound.compound_image = compound_image_path  

        try:
            db.session.commit() 
            flash('Compound updated successfully!', 'success')
            print("Update successful")
        except Exception as e:
            db.session.rollback()
            flash(f"An error occurred while updating the compound. Please try again. {e}", "danger")
            print("Error committing to database:", e)
        
        return render_template(
            'editCompound.html', 
            form=form, 
            compound=compound, 
            history_records=history_records,
            related_taxa=related_taxa, 
            logged_in=logged_in
        )

    else:
        #flash(f"An error occurred while updating the compound. Please try again. {form.errors}", "danger")
        print("Form validation errors:", form.errors)

    return render_template(
            'editCompound.html', 
            form=form, 
            compound=compound, 
            history_records=history_records,
            related_taxa=related_taxa, 
            logged_in=logged_in
        )


@compounds.route('/compound/revert/<int:history_id>', methods=['POST'])
@login_required
def revert_compound(history_id):
    history_record = CompoundHistory.query.get_or_404(history_id)
    main_compound = Compounds.query.get_or_404(history_record.compound_id)

    current_state_history = CompoundHistory(
        compound_id=main_compound.id,
        smiles=main_compound.smiles,
        journal=main_compound.journal,
        article_url=main_compound.article_url,
        inchi_key=main_compound.inchi_key,
        exact_molecular_weight=main_compound.exact_molecular_weight
    )
    db.session.add(current_state_history)

    main_compound.journal = history_record.journal
    main_compound.smiles = history_record.smiles
    main_compound.article_url = history_record.article_url
    main_compound.inchi_key = history_record.inchi_key
    main_compound.exact_molecular_weight = history_record.exact_molecular_weight
    main_compound.class_results = history_record.class_results
    main_compound.superclass_results = history_record.superclass_results
    main_compound.pathway_results = history_record.pathway_results
    main_compound.pubchem_id = history_record.pubchem_id
    main_compound.inchi = history_record.inchi

    try:
        db.session.commit()
        flash(f"Compound reverted to version from {history_record.created_at.strftime('%Y-%m-%d %H:%M')}.", "success")
    except Exception as e:
        db.session.rollback()
        flash(f"An error occurred while reverting the compound: {e}", "danger")

    return redirect(url_for('compounds.edit_compound', id=main_compound.id))

@compounds.route('/compound/<int:id>/delete', methods=['POST'])
@login_required
def delete_compound(id):
    compound = Compounds.query.get_or_404(id)

    try:
        for doi in compound.dois:
            DOI.soft_delete(doi)

        Compounds.soft_delete(compound)
        flash('Compound deleted successfully!', 'success')
        return redirect(url_for('main.index'))
    except Exception as e:
        db.session.rollback()
        flash(f"An error occurred while deleting the compound: {e}", "danger")
        return redirect(url_for('compounds.edit_compound', id=id))

@compounds.route('/compound/<int:id>/restore', methods=['POST'])
@login_required
def restore_compound(id):
    try:
        compound = Compounds.query.get_or_404(id)
        compound.restore
        return jsonify({"message": "Compound restored successfully"}), 200
    except Exception as e: 
        return jsonify({"message": "An error occurred while restoring the compound"}), 500

@compounds.route('/download_compounds', methods=['GET'])
def download_compounds():
    logged_in = current_user.is_authenticated
    compounds = Compounds.query.filter_by(status='public').all()

    if 'download' in request.args:  
        output = StringIO()
        writer = csv.DictWriter(output, fieldnames = [
        'id', 'journal', 'compound_name', 'compound_image', 'smiles', 'article_id',
        'inchi', 'inchikey', 'exactmolwt', 'pubchem', 'source', 'user_id', 
        'status', 'class_results', 'superclass_results', 'pathway_results', 'isglycoside'
        ])
        writer.writeheader()

        for compound in compounds:
            writer.writerow(compound.to_dict())

        output.seek(0)

        return Response(
            output.getvalue(),
            mimetype='text/csv',
            headers={'Content-Disposition': 'attachment; filename=compounds.csv'}
        )

    return render_template('download_compounds.html', compounds=compounds, logged_in=logged_in)

