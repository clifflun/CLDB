import pandas as pd
import numpy as np
import sqlite3 as sq
from sqlalchemy import create_engine
import glob
import os

pd.set_option('display.expand_frame_repr', False)

# --- CONFIGURATION ---
META_FILE = 'Z:/Members/clun/CLDB/meta/meta_SR_hg38.tsv' 
CNV_FILE  = 'Z:/Members/clun/CLDB/CLDB_CNV_hg38_updated_annotated.tsv'
GAP_FILE  = 'Z:/Members/clun/CLDB/data/3.0.0/CLDB_CNV_hg38_gap_diff.tsv'
DB_NAME   = 'prod_CLDB_SR.sqlite'

def write_meta_to_DB(fn):
    print(f"Loading metadata from {fn}...")
    try:
        # Added encoding='utf-8' for safety
        df = pd.read_csv(fn, sep='\t', encoding='utf-8')
        df = df[['pt_id', 'sex', 'phenotype']]
        engine = create_engine(f'sqlite:///{DB_NAME}')
        df.to_sql('meta_hg38', con=engine, if_exists='replace', index=False)
        print("Metadata loaded successfully.")
    except Exception as e:
        print(f"Error loading metadata: {e}")

def write_to_DB():
    engine = create_engine(f'sqlite:///{DB_NAME}')
    fnames = glob.glob(CNV_FILE)
    
    if not fnames:
        print(f"ERROR: Could not find CNV file: {CNV_FILE}")
        return

    for fn in fnames: 
        print(f"Loading CNV Data from: {fn}")
        try:
            # Chunking with low_memory=False to avoid DtypeWarnings
            # Added encoding='utf-8'
            with pd.read_csv(fn, sep='\t', chunksize=100000, low_memory=False, encoding='utf-8') as reader:
                for idx, chunk in enumerate(reader):
                    print(f'  Processing chunk {idx}...')
                    chunk.to_sql('CNV_hg38', con=engine, if_exists='append', index=False)
            print("CNV Data loaded successfully.")
        except Exception as e:
            print(f"Error loading CNV data: {e}")

def create_index():
    print("Creating indexes for CNV_hg38...")
    conn = sq.connect(DB_NAME)
    cursor = conn.cursor()
    _list = ['chrom1', 'chrom2', 'SV_ID', 'SV_TYPE', 'PT_ID', 'FAM_ID', 'PROJECT', 'IS_PROBAND', 'CLUSTER_ID']
    
    for i in _list:
        try:
            print(f"  Creating index for column: {i}")
            cursor.execute(f'CREATE INDEX if NOT EXISTS {i}_index ON CNV_hg38 ({i})')
        except Exception as e:
            print(f"  Skipping index for {i} (column may not exist)")
            pass
            
    conn.commit()
    conn.close()
    print("Indexes created.")

def delete():
    print('Deleting old CNV_hg38 table...')
    conn = sq.connect(DB_NAME)
    cursor = conn.cursor()    
    cursor.execute('DROP TABLE IF EXISTS CNV_hg38')
    conn.commit()
    conn.close()
    print('Table deleted.')

def update_column_optimized():
    print(f"Starting optimized update for SD_overlap in CNV...")
    
    if not os.path.exists(GAP_FILE):
        print(f"WARNING: Gap file not found at {GAP_FILE}. Skipping update.")
        return

    try:
        # 1. Load gap data
        print("  Loading gap file...")
        df = pd.read_csv(GAP_FILE, sep='\t', encoding='utf-8')
        conn = sq.connect(DB_NAME)
        
        # 2. Load into temp table
        print("  Writing to temporary table 'gap_update'...")
        df.to_sql('gap_update', conn, if_exists='replace', index=False)
        
        cursor = conn.cursor()
        
        # 3. Create Index on Temp Table (Critical for Speed)
        print("  Indexing temp table...")
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_gap_uuid ON gap_update(UUID)')
        
        # 4. Create Index on Target Table (Critical for Speed)
        # Note: Previous checks confirmed CNV_hg38 has a UUID column.
        print("  Indexing target table (UUID)...")
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_cnv_uuid ON CNV_hg38(UUID)') 
        
        # 5. Execute Update
        # Corrected SQL: References CNV_hg38 correctly
        print("  Running SQL Update...")
        update_query = """
        UPDATE CNV_hg38
        SET SD_overlap = (
            SELECT SD_overlap
            FROM gap_update
            WHERE gap_update.UUID = CNV_hg38.UUID
        )
        WHERE UUID IN (SELECT UUID FROM gap_update);
        """
        cursor.execute(update_query)
        conn.commit()
        print(f"  Update complete. Rows affected: {cursor.rowcount}")
        
    except Exception as e:
        print(f"ERROR during update: {e}")
    finally:
        if conn: conn.close()

def query():
    print("Verifying database content...")
    conn = sq.connect(DB_NAME)
    try:
        df = pd.read_sql("SELECT COUNT(*) as Count FROM CNV_hg38", conn)
        count = df.iloc[0]['Count']
        print(f"Total Rows in CNV_hg38: {count}")
        if count == 0:
             print("WARNING: Table is empty.")
    except Exception as e:
        print(f"Query failed: {e}")
    conn.close()

def main():
    print("--- Starting CNV Database Update ---")
    delete()
    write_to_DB()
    create_index()
    write_meta_to_DB(META_FILE)
    update_column_optimized() 
    query()
    print("--- CNV Update Complete ---")

if __name__ == '__main__':
    main()