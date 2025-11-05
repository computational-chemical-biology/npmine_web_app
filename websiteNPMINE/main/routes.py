from flask import Blueprint, render_template, request, jsonify, abort, redirect, url_for,flash
from flask_login import login_required, current_user, login_user, logout_user
from websiteNPMINE.models import Compounds, DOI, Accounts, Taxa, doicomp, doitaxa
from websiteNPMINE import db
import requests
from sqlalchemy import or_
from sqlalchemy.orm.attributes import InstrumentedAttribute
from sqlalchemy.orm import aliased
import collections
from websiteNPMINE import csrf
from rdkit import Chem
from sqlalchemy.sql import func
from sqlalchemy import or_, desc, asc, text

main = Blueprint('main', __name__)

@main.route('/')
@main.route("/home")
def home():
    logged_in = current_user.is_authenticated  
    return render_template("index.html", logged_in=logged_in)

@main.route('/api/benchmark')
def benchmark_data():
    """
    Executa um benchmark de performance comparando a busca Tanimoto
    lenta (full scan) contra a rápida (índice BFP).
    """
    smiles = request.args.get('smiles')
    
    if not smiles:
        return jsonify({"erro": "Parâmetro 'smiles' é obrigatório. (ex: ?smiles=C[C@H]1...)"}), 400
    
    query_mol = Chem.MolFromSmiles(smiles)
    if not (query_mol and query_mol.GetNumHeavyAtoms() > 1):
        return jsonify({"erro": "SMILES inválido ou simples demais para o benchmark"}), 400

    slow_sql = text("""
        EXPLAIN (ANALYZE, FORMAT JSON)
        SELECT COUNT(*)
        FROM compounds
        WHERE deleted_at IS NULL
          AND tanimoto_sml(
                morganbv_fp(molecule),
                morganbv_fp(mol_from_smiles(:smiles_param))
              ) > 0.7;
    """)

    fast_sql = text("""
        EXPLAIN (ANALYZE, FORMAT JSON)
        SELECT COUNT(*)
        FROM compounds
        WHERE deleted_at IS NULL
          AND fingerprint % morganbv_fp(mol_from_smiles(:smiles_param));
    """)

    try:
        slow_result = db.session.execute(slow_sql, {"smiles_param": smiles}).fetchone()
        
        slow_plan = slow_result[0][0]['Plan']
        slow_time_ms = slow_plan['Actual Total Time']

        db.session.execute(text("SET rdkit.tanimoto_threshold = 0.7"))
        
        fast_result = db.session.execute(fast_sql, {"smiles_param": smiles}).fetchone()
        
        fast_plan = fast_result[0][0]['Plan']
        fast_time_ms = fast_plan['Actual Total Time']
        
        return jsonify({
            "status": "Benchmark Concluído",
            "composto_teste": smiles,
            "metodo_lento_ms": slow_time_ms,
            "metodo_rapido_ms": fast_time_ms,
            "melhoria_de_performance": f"{(slow_time_ms / fast_time_ms):.2f}x mais rápido",
            "plano_execucao_lento": slow_plan,
            "plano_execucao_rapido": fast_plan
        })

    except Exception as e:
        db.session.rollback() 
        return jsonify({"erro": f"Falha no benchmark: {str(e)}"}), 500
    
@main.route('/api/data')
def data():
    start = request.args.get('start', type=int, default=0)
    length = request.args.get('length', type=int, default=10)
    show_deleted = request.args.get('show_deleted', 'false').lower() == 'true'
    search_term = request.args.get('search', '').strip()
    sort_str = request.args.get('sort', '')
    is_glycoside_filter = request.args.get('is_glycoside', 'all').lower()

    query_mol = Chem.MolFromSmiles(search_term) if search_term else None
    is_similarity_search = query_mol and query_mol.GetNumHeavyAtoms() > 1

    base_query = db.session.query(Compounds)

    delete_filter = Compounds.deleted_at.is_not(None) if show_deleted else Compounds.deleted_at.is_(None)
    base_query = base_query.filter(delete_filter)

    if current_user.is_authenticated:
        ownership_filter = or_(
            Compounds.ispublic == True,
            Compounds.user_id == current_user.id
        )
        base_query = base_query.filter(ownership_filter)
    else:
        base_query = base_query.filter(Compounds.ispublic == True)

    if is_glycoside_filter == 'yes':
        base_query = base_query.filter(Compounds.isglycoside == True)
    elif is_glycoside_filter == 'no':
        base_query = base_query.filter(Compounds.isglycoside == False)

    order = []
    similarity_col = None

    if is_similarity_search:
        
        similarity_threshold = 0.7  
        db.session.execute(text(f"SET rdkit.tanimoto_threshold = {similarity_threshold}"))
        
        query_mol_sql = func.mol_from_smiles(search_term)
        fp_query = func.morganbv_fp(query_mol_sql) 

        similarity_col = func.tanimoto_sml(Compounds.fingerprint, fp_query).label('similarity')

        search_filter = Compounds.fingerprint.op('%')(fp_query)
                
        base_query = base_query.filter(search_filter)
        order.append(desc('similarity'))

    elif search_term:   
        print(f"Running Text Search for: {search_term}")
        
        search_filter_text = f"%{search_term}%"
        search_filter = or_(
            Compounds.compound_name.ilike(search_filter_text),
            Compounds.inchi_key.ilike(search_filter_text),
            Compounds.pubchem_id.ilike(search_filter_text)
        )
        base_query = base_query.filter(search_filter)
        
        sort_columns = ['id', 'compound_name', 'inchi_key', 'pubchem_id']
        if sort_str:
            for field in sort_str.split(','):
                if not field: continue
                direction = 'desc' if field.startswith('-') else 'asc'
                name = field.lstrip('-')
                
                if name in sort_columns:
                    col = getattr(Compounds, name)
                    order.append(getattr(col, direction)())

    if not order: 
        order.append(Compounds.id.asc())

    count_query = base_query.with_entities(func.count(Compounds.id))
    total = count_query.scalar()

    data_query = None
    if is_similarity_search:
        data_query = base_query.with_entities(Compounds, similarity_col)
    else:
        data_query = base_query.with_entities(Compounds)
    
    data_query = data_query.order_by(*order).offset(start).limit(length)
    results = data_query.all()

    data = []
    for row in results:
        compound = None
        compound_data = {}

        if is_similarity_search:
            compound = row.Compounds
            compound_data = compound.to_dict()
            compound_data['similarity'] = round(row.similarity, 3) 
        else:
            compound = row
            compound_data = compound.to_dict()

        compound_data['compound_image'] = url_for('static', filename=compound.compound_image) if compound.compound_image else None
        data.append(compound_data)

    return jsonify({'data': data, 'total': total})


@main.route('/compound/<int:compound_id>')
def compound(compound_id):
    compound = Compounds.query.get(compound_id)
    logged_in = current_user.is_authenticated  
    if not compound:
        abort(404)

    articles = [(doi.id, doi.doi) for doi in compound.dois]

    # Use NPclassifier API
    #api_url = f'https://npclassifier.gnps2.org/classify?smiles={compound.smiles}'
    #resposta_api = requests.get(api_url).json()

    return render_template('compound.html', logged_in=logged_in, compound=compound, articles=articles)


@main.route('/article/<int:article_id>')
def article(article_id):
    logged_in = current_user.is_authenticated  
    
    doi_record = DOI.query.get(article_id)
    if not doi_record:
        abort(404)
    
    if logged_in:
        compounds = [compound for compound in doi_record.compounds 
                     if compound.status == 'public' or compound.user_id == current_user.id]
    else:
        compounds = [compound for compound in doi_record.compounds 
                     if compound.status == 'public']

    if compounds:
        first_compound = compounds[0]
        journal_name = first_compound.journal
        article_url = first_compound.article_url
        created_at = first_compound.created_at
    else:
        journal_name = None
        article_url = None
        created_at = None

    related_compound_names = []
    compound_name_to_id_map = {}
    verbatim_values = set()

    taxa_for_doi = Taxa.query.filter(Taxa.dois.any(id=article_id)).all()
    taxa_verbatim_values = set(taxon.verbatim for taxon in taxa_for_doi)

    for compound in compounds:
        compound_name = compound.compound_name
        compound_id = compound.id
        related_compound_names.append(compound_name)
        compound_name_to_id_map[compound_name] = compound_id

        for taxon in taxa_for_doi:
            if taxon in doi_record.taxa:
                verbatim_values.add(taxon.verbatim) 

    return render_template('article.html', 
                           doi_record=doi_record, 
                           journal_name=journal_name, 
                           article_url=article_url, 
                           created_at=created_at, 
                           related_compound_names=related_compound_names, 
                           compound_name_to_id_map=compound_name_to_id_map, 
                           verbatim_values=verbatim_values, 
                           logged_in=logged_in)



@main.route('/profile/<int:profile_id>')
@login_required
def profile(profile_id):
    logged_in = current_user.is_authenticated
    user = Accounts.query.get(profile_id)
    if not user:
        abort(404) 

    public_compounds = Compounds.query.filter_by(user_id=profile_id, status='public').all()
    private_compounds = Compounds.query.filter_by(user_id=profile_id, status='private').all()

    return render_template('profile.html', logged_in=logged_in, user=user, public_compounds=public_compounds, private_compounds=private_compounds)

def delete_compound_and_related(compound):
    for doi in compound.dois:
        if len(doi.compounds) == 1:
            for taxon in doi.taxa:
                taxon_related_compounds = sum(
                    len(d.compounds) for d in taxon.dois if d != doi
                )
                if taxon_related_compounds == 0:
                    db.session.delete(taxon)

            DOI.soft_delete(doi)

    Compounds.soft_delete(compound)

@main.route('/compound/<int:compound_id>/delete', methods=['POST'])
@csrf.exempt
@login_required
def delete_compound(compound_id):
    compound = Compounds.query.get_or_404(compound_id)

    if not current_user.role_id == 1 and current_user.id != compound.user_id:
        flash("You don't have permission to delete this compound.", "error")
        return redirect(url_for('main.home'))

    try:
        delete_compound_and_related(compound)
        db.session.commit()
        flash("Compound deleted successfully.", "success")
    except Exception as e:
        db.session.rollback()
        flash(f"An error occurred while deleting the compound: {str(e)}", "error")

    return redirect(url_for('main.home'))


@main.route('/toggle_privacy/<int:compound_id>', methods=['POST'])
@csrf.exempt
@login_required
def toggle_privacy(compound_id):
    compound = Compounds.query.get_or_404(compound_id)

    if current_user.id != compound.user_id:
        abort(403)  
    if compound.status == 'private':
        compound.status = 'public'
        flash(f"Compound '{compound.compound_name}' is now public.", "success")
    else:
        compound.status = 'private'
        flash(f"Compound '{compound.compound_name}' is now private.", "success")

    db.session.commit()


    return redirect(url_for('main.profile', profile_id=current_user.id))




