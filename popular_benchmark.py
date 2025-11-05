import pandas as pd
import psycopg2
from tqdm import tqdm
from datetime import datetime
import os
from dotenv import load_dotenv

load_dotenv()

DB_HOST = os.getenv('DB_HOST')
DB_PORT = os.getenv('DB_PORT')
DB_NAME = os.getenv('DB_NAME')
DB_USER = os.getenv('DB_USER')
DB_PASSWORD = os.getenv('DB_PASSWORD')
NPMINE_WEB_APP_PASSWORD = os.getenv('NPMINE_WEB_APP_PASSWORD')
NPMINE_WEB_APP_EMAIL = os.getenv('NPMINE_WEB_APP_EMAIL')

ARQUIVO_DADOS = 'coconut_csv_lite-11-2025.csv' 
NUMERO_COMPOSTOS = 10000 

COLUNAS_CSV_PARA_LER = [
    'name',
    'canonical_smiles',     
    'standard_inchi_key', 
    'exact_molecular_weight',
    'standard_inchi',       
]

def conectar_banco():
    """Conecta ao banco de dados PostgreSQL."""
    try:
        conn = psycopg2.connect(
            host=DB_HOST,
            port=DB_PORT,
            database=DB_NAME,
            user=DB_USER,
            password=DB_PASSWORD
        )
        print("Conexão ao PostgreSQL bem-sucedida.")
        return conn
    except psycopg2.OperationalError as e:
        print(f"Erro ao conectar ao banco: {e}")
        return None

def get_admin_id(cursor):
    """Pega o ID do primeiro usuário admin (role_id = 1)."""
    try:
        cursor.execute("SELECT id FROM accounts WHERE role_id = 1 LIMIT 1")
        admin = cursor.fetchone()
        if admin:
            print(f"Admin ID encontrado: {admin[0]}")
            return admin[0]
        else:
            print("Nenhum usuário Admin (role_id=1) encontrado. Usando ID 1 como padrão.")
            return 1 
    except Exception as e:
        print(f"Erro ao buscar admin ID: {e}")
        return 1

def popular_dados():
    conn = conectar_banco()
    if not conn:
        return
    
    cursor = conn.cursor()
    admin_id = get_admin_id(cursor)

    try:
        print(f"Lendo o arquivo {ARQUIVO_DADOS}...")
        df = pd.read_csv(
            ARQUIVO_DADOS, 
            usecols=COLUNAS_CSV_PARA_LER,
            nrows=NUMERO_COMPOSTOS,
            low_memory=False 
        )
        df = df.dropna(subset=['canonical_smiles', 'standard_inchi_key'])
        print(f"Encontrados {len(df)} compostos válidos para inserir.")

    except FileNotFoundError:
        print(f"Erro: Arquivo {ARQUIVO_DADOS} não encontrado.")
        return
    except KeyError as e:
        print(f"Erro: Coluna não encontrada no CSV: {e}")
        return
    except Exception as e:
        print(f"Erro ao ler o CSV: {e}")
        return

    for index, row in tqdm(df.iterrows(), total=df.shape[0], desc="Populando Banco"):
        
        compound_name_raw = row.get('name')
        compound_name = None if pd.isna(compound_name_raw) else compound_name_raw
        smiles = row.get('canonical_smiles')
        inchi_key = row.get('standard_inchi_key')
        inchi = row.get('standard_inchi')
        
        mw = row.get('exact_molecular_weight')
        exact_molecular_weight = float(mw) if pd.notna(mw) and isinstance(mw, (int, float, str)) and str(mw).replace('.', '', 1).isdigit() else None

        article_url = None
        pubchem_id = None
        class_results = None
        superclass_results = None
        pathway_results = None
        isglycoside = None

        journal = None
        compound_image = None
        ispublic = False
        source = 'COCONUT'
        status = 'private'
        created_at = datetime.utcnow()

        try:
            cursor.execute(
                """
                INSERT INTO compounds (
                    journal, compound_name, compound_image, smiles, article_url, 
                    inchi_key, exact_molecular_weight, class_results, superclass_results, 
                    pathway_results, isglycoside, ispublic, pubchem_id, inchi, 
                    source, user_id, status, created_at
                ) 
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                """,
                (
                    journal, compound_name, compound_image, smiles, article_url,
                    inchi_key, exact_molecular_weight, class_results, superclass_results,
                    pathway_results, isglycoside, ispublic, pubchem_id, inchi,
                    source, admin_id, status, created_at
                )
            )
        except Exception as e:
            print(f"Erro ao inserir InChIKey {inchi_key}: {e}")
            conn.rollback() 
            continue
    
    try:
        conn.commit()
        print(f"\nSucesso! {len(df)} compostos processados e inseridos.")
    except Exception as e:
        print(f"Erro ao commitar as alterações: {e}")
        conn.rollback()
    finally:
        cursor.close()
        conn.close()

if __name__ == "__main__":
    popular_dados()