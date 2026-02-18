from . import db, login_manager
from flask_login import UserMixin
from itsdangerous import URLSafeTimedSerializer as Serializer
from flask import current_app
from datetime import datetime
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy.types import UserDefinedType


admin_role_id = 1
editor_role_id = 2
user_role_id = 3

@login_manager.user_loader
def load_user(user_id):
    return Accounts.query.get(int(user_id))

class MOL(UserDefinedType):
    def get_col_spec(self, **kw):
        return "MOL"
    
class BFP(UserDefinedType):
    """
    Custom SQLAlchemy type for the RDKit 'BFP' (bit fingerprint) data type.
    """
    def get_col_spec(self, **kw):
        return "BFP"

class Role(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(50), unique=True, nullable=False)

    def __repr__(self):
        return f'<Role: {self.name}>'

class Accounts(db.Model, UserMixin):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(20), unique=True, nullable=False)
    name = db.Column(db.String(256), nullable=False)
    surname = db.Column(db.String(256), nullable=False)
    academic_level = db.Column(db.String(32), nullable=True)
    email = db.Column(db.String(120), unique=True, nullable=False)
    password = db.Column(db.String(64), nullable=False)
    created_at = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
    updated_at = db.Column(db.DateTime, nullable=False, default=datetime.utcnow, onupdate=datetime.utcnow)
    deleted_at = db.Column(db.DateTime, nullable=True)

    role_id = db.Column(db.Integer, db.ForeignKey('role.id'), nullable=False, default=user_role_id)  
    role = db.relationship('Role', backref=db.backref('accounts', lazy='dynamic'))
    
    def get_reset_token(self, expires_sec=1800):
        s = Serializer(current_app.config['SECRET_KEY'])
        return s.dumps({'user_id': self.id})

    @staticmethod
    def verify_reset_token(token):
        s = Serializer(current_app.config['SECRET_KEY'])
        try:
            user_id = s.loads(token)['user_id']
        except:
            return None
        return Accounts.query.get(user_id)

    def __repr__(self):
        return f"User('{self.username}', '{self.email}')"
    
    def soft_delete(self):
        self.deleted_at = datetime.utcnow()

    def restore(self):
        self.deleted_at = None

# Association tables with delete cascade
doicomp = db.Table(
    'doicomp',
    db.Column('doi_id', db.Integer, db.ForeignKey('doi.id', ondelete="CASCADE"), primary_key=True),
    db.Column('compound_id', db.Integer, db.ForeignKey('compounds.id', ondelete="CASCADE"), primary_key=True)
)

doitaxa = db.Table(
    'doitaxa',
    db.Column('doi_id', db.Integer, db.ForeignKey('doi.id', ondelete="CASCADE"), primary_key=True),
    db.Column('taxon_id', db.Integer, db.ForeignKey('taxa.id', ondelete="CASCADE"), primary_key=True)
)

class DOI(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    doi = db.Column(db.String(150), unique=True)

    compounds = db.relationship('Compounds', 
                                secondary=doicomp, 
                                back_populates='dois')
    taxa = db.relationship('Taxa', 
                           secondary=doitaxa, 
                           back_populates='dois')

    def __repr__(self):
        return f'<DOI: {self.doi}>'
    
    def soft_delete(self):
        self.deleted_at = datetime.utcnow()
        db.session.add(self)
        db.session.commit()

    def restore(self):
        self.deleted_at = None
        db.session.add(self)
        db.session.commit()

class Compounds(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    journal = db.Column(db.String(5000))
    compound_name = db.Column(db.String(5000))
    compound_image = db.Column(db.String(5000))
    molecule = db.Column(MOL(), nullable=True)
    smiles = db.Column(db.String(5000))
    fingerprint = db.Column(BFP(), nullable=True)
    article_url = db.Column(db.String(500))
    inchi_key = db.Column(db.String(5000))
    exact_molecular_weight = db.Column(db.Float)
    class_results = db.Column(db.String(5000))
    superclass_results = db.Column(db.String(5000))
    pathway_results = db.Column(db.String(5000))
    isglycoside = db.Column(db.Boolean, default=True, nullable=True)
    ispublic = db.Column(db.Boolean, default=False, nullable=True)
    pubchem_id = db.Column(db.String(5000))
    inchi = db.Column(db.String(5000))
    source = db.Column(db.String(10))
    user_id = db.Column(db.Integer, db.ForeignKey('accounts.id'))
    status = db.Column(db.String(10), nullable=False, default='private')
    created_at = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
    updated_at = db.Column(db.DateTime, nullable=True, default=datetime.utcnow)
    deleted_at = db.Column(db.DateTime, nullable=True)

    account = db.relationship("Accounts", backref="compounds")
    dois = db.relationship('DOI', 
                           secondary=doicomp, 
                           back_populates='compounds')

    def __repr__(self):
        return f'Compounds: {self.id}'

    def to_dict(self):
        return {
            'id': self.id,
            'journal': self.journal,
            'compound_name': self.compound_name,
            'compound_image': self.compound_image,
            'smiles': self.smiles,
            'article_id': self.article_url,
            'inchi': self.inchi,
            'inchikey': self.inchi_key,
            'exactmolwt': self.exact_molecular_weight,
            'pubchem': self.pubchem_id,
            'source': self.source,
            'user_id': self.user_id,
            'status': self.status,
            'class_results': self.class_results,
            'superclass_results': self.superclass_results,
            'pathway_results': self.pathway_results,
            'isglycoside': self.isglycoside,
        }
    
    def soft_delete(self):
        self.deleted_at = datetime.utcnow()
        db.session.add(self)
        db.session.commit()

    def restore(self):
        self.deleted_at = None
        db.session.add(self)
        db.session.commit()

    history = db.relationship('CompoundHistory', 
        backref='compound', 
        lazy='dynamic', 
        cascade='all, delete-orphan')

class Taxa(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    article_url = db.Column(db.String(5000))
    verbatim = db.Column(db.String(5000))
    odds = db.Column(db.Float)
    datasourceid = db.Column(db.String(5000))
    taxonid = db.Column(db.Integer)
    classificationpath = db.Column(db.String(5000))
    classificationrank = db.Column(db.String(5000))
    matchtype = db.Column(db.String(50))
    user_id = db.Column(db.Integer, db.ForeignKey('accounts.id'))
    created_at = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
    
    account = db.relationship("Accounts", backref="taxons")
    dois = db.relationship('DOI', 
                           secondary=doitaxa, 
                           back_populates='taxa')

    def __repr__(self):
        return f'<Taxa: {self.verbatim}, {self.article_url}, {self.classificationpath}, {self.classificationrank}, {self.user_id}>'
    
class CompoundHistory(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    compound_id = db.Column('compound_id', db.Integer, db.ForeignKey('compounds.id', ondelete="CASCADE"))
    active = db.Column(db.Boolean, default=True, nullable=True) 
    journal = db.Column(db.String(5000))
    smiles = db.Column(db.String(5000))
    article_url = db.Column(db.String(500))
    inchi_key = db.Column(db.String(5000))
    exact_molecular_weight = db.Column(db.Float)
    class_results = db.Column(db.String(5000))
    superclass_results = db.Column(db.String(5000))
    pathway_results = db.Column(db.String(5000))
    pubchem_id = db.Column(db.String(5000))
    inchi = db.Column(db.String(5000))
    source = db.Column(db.String(10))
    created_at = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)