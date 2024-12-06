from flask import Blueprint, render_template, flash, redirect, url_for, request, current_app, g,jsonify, Response
from flask_login import login_required, current_user
from websiteNPMINE.compounds.forms import CompoundForm, SearchForm, CompoundEditForm
from websiteNPMINE.models import Compounds,DOI,Taxa
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
    # Generate the filename
    filename = f"{compound_id}.png"
    relative_path = os.path.join('compound_images', filename)
    
    # Normalize the path and replace backslashes with forward slashes
    relative_path = os.path.normpath(relative_path).replace("\\", "/")

    # Define the full path to save the image
    full_path = os.path.join(current_app.root_path, 'static', relative_path)

    # Ensure the directory exists
    os.makedirs(os.path.dirname(full_path), exist_ok=True)

    # Check if the image already exists; if not, generate it
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol, size=(200, 200))
    img.save(full_path)

    # Return the relative path to be stored in the database
    return relative_path

@compounds.route('/new_compound', methods=['GET', 'POST'])
@login_required
def registerCompound():
    logged_in = current_user.is_authenticated
    form = CompoundForm()

    if form.validate_on_submit():
        doi = form.doi.data

        # Ensure DOI is provided
        if not doi:
            flash('DOI Link is required', 'error')
            return redirect(request.url)

        # Check if DOI exists in the DOI table
        existing_doi = DOI.query.filter_by(doi=doi).first()
        if not existing_doi:
            new_doi = DOI(doi=doi)
            db.session.add(new_doi)
            db.session.commit()
            existing_doi = new_doi
        else:
            flash('DOI already in database!', 'info')

        compound_blocks = request.form.getlist('inchikey')
        for i, inchikey in enumerate(compound_blocks):
            genus = request.form.getlist('genus')[i]
            species = request.form.getlist('species')[i]
            origin_type = request.form.getlist('origin_type')[i]  # Get the origin type

            if not inchikey:
                flash(f'InChI Key is required for Compound {i + 1}', 'error')
                continue

            # Check if the compound already exists by InChI Key
            existing_compound = Compounds.query.filter_by(inchi_key=inchikey).first()

            if existing_compound:
                # Associate the existing compound with the new DOI if not already associated
                if existing_doi not in existing_compound.dois:
                    existing_compound.dois.append(existing_doi)
                    db.session.commit()
                flash(f'Compound with InChI Key {inchikey} already in database; associated with new DOI.', 'info')
            else:
                pubchem_data = fetch_pubchem_data(inchikey)
                if not pubchem_data:
                    flash(f'Failed to fetch data from PubChem for Compound {i + 1}', 'error')
                    continue

                # Create a new compound
                compound = Compounds(
                    journal=None,
                    compound_name=pubchem_data['compound_name'],
                    smiles=pubchem_data['smiles'],
                    article_url=doi,
                    inchi_key=inchikey,
                    exact_molecular_weight=pubchem_data['exact_molecular_weight'],
                    class_results=None,
                    superclass_results=None,
                    pathway_results=None,
                    isglycoside=None,
                    pubchem_id=pubchem_data['pubchem_id'],
                    inchi=pubchem_data['inchi'],
                    source='NPMine',
                    user_id=current_user.id
                )
                db.session.add(compound)
                db.session.commit()

                # Associate the compound with the DOI
                compound.dois.append(existing_doi)
                db.session.commit()

                # Save compound image
                smiles = pubchem_data.get('smiles')
                img_path = save_compound_image(compound.id, smiles)
                if img_path:
                    compound.compound_image = img_path
                    db.session.commit()

            if genus and species:
                # Check if the taxa already exists
                existing_taxon = Taxa.query.filter_by(verbatim=f"{genus} {species}").first()
                if not existing_taxon:
                    # Create new taxa if it does not exist
                    taxon = Taxa(
                        article_url=doi,
                        verbatim=f"{genus} {species}",
                        user_id=current_user.id
                    )
                    db.session.add(taxon)
                    db.session.commit()

                    # Associate the new taxon with the DOI
                    taxon.dois.append(existing_doi)
                else:
                    # Associate the existing taxon with the DOI if not already associated
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
        # Search for the compound using the InChI Key
        compound = pcp.get_compounds(inchikey, 'inchikey')[0]
        return {
            'compound_name': compound.synonyms[0] if compound.synonyms else None,
            'smiles': compound.canonical_smiles,
            'inchi': compound.inchi,
            'exact_molecular_weight': compound.molecular_weight,
            'pubchem_id': compound.cid
        }
    except pcp.PubChemHTTPError as e:
        print(f"PubChem HTTP Error: {e}")
        return {}
    except pcp.PubChemPyError as e:
        print(f"PubChemPy Error: {e}")
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
                (Compounds.status == 'public')  # Only show Public compounds
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
        # Perform a query on the DOI model, joining with the Compounds model
        results = (
            DOI.query
            .join(DOI.compounds)  # Join DOI with Compounds through the relationship
            .filter(Compounds.status == 'public')  # Filter for public compounds only
            .filter(DOI.doi.ilike(f"%{q}%"))
            .distinct()
            .options(
                db.joinedload(DOI.compounds),
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
        # Query for Taxa associated with DOIs that have public compounds
        results = (
            Taxa.query
            .join(Taxa.dois)  # Join Taxa with DOI through the relationship
            .join(DOI.compounds)  # Join DOI with Compounds to access status
            .filter(Compounds.status == 'public')  # Filter to only public compounds
            .filter(Taxa.verbatim.ilike(f"%{q}%"))
            .distinct()
            .options(
                db.joinedload(Taxa.dois).joinedload(DOI.compounds),
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

# Function to process similarity search
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
    
    # Sort results by similarity score in descending order
    search_res.sort(key=lambda x: x[2], reverse=True)
    
    return search_res

@compounds.route('/search_structure', methods=['GET', 'POST'])
@csrf.exempt
def search_structure():
    logged_in = current_user.is_authenticated

    if request.method == 'GET':
        return render_template('search_structure.html')

    # Retrieve only public compounds from the database
    cpds = db.session.query(Compounds).filter(Compounds.status == 'public').all()
    search_res = []
    query = None

    # Handle text SMILES input
    if 'textSmiles' in request.form and request.form['textSmiles']:
        smiles = request.form['textSmiles']
        try:
            inchi = Chem.MolToInchi(Chem.MolFromSmiles(smiles))
            query = f"https://cactus.nci.nih.gov/chemical/structure/{quote(inchi)}/image"
            print(f"Generated Image URL: {query}")  # Log the generated URL for debugging
            search_res = process_similarity_search(inchi, cpds)
        except Exception as e:
            print(f"Error processing SMILES input: {e}")
            return "Invalid SMILES input", 400

    # Handle drawn SMILES input
    elif 'drawSmiles' in request.form and request.form['drawSmiles']:
        smiles = request.form['drawSmiles']
        try:
            inchi = Chem.MolToInchi(Chem.MolFromSmiles(smiles))
            query = f"https://cactus.nci.nih.gov/chemical/structure/{quote(inchi)}/image"
            print(f"Generated Image URL: {query}")  # Log the generated URL for debugging
            search_res = process_similarity_search(inchi, cpds)
        except Exception as e:
            print(f"Error processing drawn SMILES: {e}")
            return "Invalid SMILES input", 400

    # Handle file upload for SMILES
    elif 'data' in request.files:
        data = request.files['data']
        if data.filename == '':
            return "No selected file", 400

        # Process uploaded file content
        try:
            file_content = data.read().decode('utf-8')
            smiles_list = file_content.splitlines()
        except Exception as e:
            print(f"Error reading file: {e}")
            return "Error reading file", 400

        # Process each SMILES in the file
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

    # No valid input found
    else:
        return "Invalid request", 400

    # Directly pass search_res to the template without using json.dumps()
    return render_template('search_results_structure.html', dfinal=search_res, query=query, logged_in=logged_in)


def update_compound_relationships(compound, doi_data_list):
    # Iterate through existing DOIs in the compound
    for existing_doi in compound.dois:
        # Find matching DOI data in the new data list
        matching_doi_data = next((doi_data for doi_data in doi_data_list if doi_data['doi'] == existing_doi.doi), None)
        if matching_doi_data:
            # Update existing DOI with new data
            existing_doi.some_attribute = matching_doi_data.get('some_attribute', existing_doi.some_attribute)
            # Remove updated DOI data from the list to handle additions separately
            doi_data_list.remove(matching_doi_data)
    
    # Now, add any remaining DOIs that are new
    for new_doi_data in doi_data_list:
        new_doi = DOI(**new_doi_data)
        compound.dois.append(new_doi)
    
    db.session.commit()

@compounds.route('/compound/<int:id>/edit', methods=['GET', 'POST'])
@login_required
def edit_compound(id):
    logged_in = current_user.is_authenticated
    compound = Compounds.query.get_or_404(id)
    form = CompoundEditForm(obj=compound)

    # Clear and repopulate DOI form fields
    form.dois.entries.clear()
    for doi in compound.dois:
        form.dois.append_entry({'doi': doi.doi})

    # Fetch related taxa
    related_taxa = [taxa for doi in compound.dois for taxa in doi.taxa]

    if form.validate_on_submit():
        # Log form data to check what is being submitted
        print("Form data:", form.data)

        # Store old SMILES for comparison
        old_smiles = compound.smiles

        # Update compound fields with form data
        compound.journal = form.journal.data or compound.journal
        compound.compound_name = form.compound_name.data or compound.compound_name
        compound.smiles = form.smiles.data or compound.smiles
        compound.exact_molecular_weight = form.exact_molecular_weight.data or compound.exact_molecular_weight
        compound.class_results = form.class_results.data or compound.class_results
        compound.superclass_results = form.superclass_results.data or compound.superclass_results
        compound.pathway_results = form.pathway_results.data or compound.pathway_results
        compound.isglycoside = form.isglycoside.data or compound.isglycoside
        compound.pubchem_id = form.pubchem_id.data or compound.pubchem_id
        compound.inchi = form.inchi.data or compound.inchi
        compound.article_url = form.article_url.data or compound.article_url
        compound.status = form.status.data or compound.status

        # Debug: Show which DOIs are being updated or added
        print("DOIs before update:", [d.doi for d in compound.dois])
        
        form_doi_data_list = [entry.data['doi'] for entry in form.dois.entries]

        # Update existing DOIs and remove those not in the form
        for existing_doi in compound.dois[:]:
            if existing_doi.doi in form_doi_data_list:
                form_doi_entry = next(entry for entry in form.dois.entries if entry.data['doi'] == existing_doi.doi)
                # Update attributes as needed
            else:
                compound.dois.remove(existing_doi)
                db.session.delete(existing_doi)

        # Add new DOIs from form
        for form_doi_entry in form.dois.entries:
            form_doi_data = form_doi_entry.data
            if not any(doi.doi == form_doi_data['doi'] for doi in compound.dois):
                new_doi = DOI(doi=form_doi_data['doi'])
                compound.dois.append(new_doi)
                db.session.add(new_doi)

        # Check if SMILES changed and regenerate image if it did
        if compound.smiles != old_smiles:
            print("SMILES changed, updating image...")
            # Regenerate the compound image
            compound_image_path = save_compound_image(compound.id, compound.smiles)
            compound.compound_image = compound_image_path  # Assuming the image path is stored in `compound_image` field

        # Commit changes to the database
        try:
            db.session.commit()  # No need to add compound again; it's already being tracked
            flash('Compound updated successfully!', 'success')
            print("Update successful")
        except Exception as e:
            db.session.rollback()
            flash("An error occurred while updating the compound. Please try again.", "danger")
            print("Error committing to database:", e)
        
        return redirect(url_for('main.compound', compound_id=compound.id))

    # If form didn't validate, print errors
    else:
        print("Form validation errors:", form.errors)

    return render_template('editCompound.html', form=form, compound=compound, related_taxa=related_taxa, logged_in=logged_in)

@compounds.route('/download_compounds', methods=['GET'])
def download_compounds():
    # Query all publicly available compounds
    compounds = Compounds.query.filter_by(status='public').all()

    if 'download' in request.args:  # Check if the user triggered the download
        # Create a CSV in memory
        output = StringIO()
        writer = csv.DictWriter(output, fieldnames = [
        'id', 'journal', 'compound_name', 'compound_image', 'smiles', 'article_id',
        'inchi', 'inchikey', 'exactmolwt', 'pubchem', 'source', 'user_id', 
        'status', 'class_results', 'superclass_results', 'pathway_results', 'isglycoside'
        ])
        writer.writeheader()

        # Write compound data to the CSV
        for compound in compounds:
            writer.writerow(compound.to_dict())

        output.seek(0)

        # Create a Response to send the CSV file as a downloadable attachment
        return Response(
            output.getvalue(),
            mimetype='text/csv',
            headers={'Content-Disposition': 'attachment; filename=compounds.csv'}
        )

    # Render the template when not downloading
    return render_template('download_compounds.html', compounds=compounds)