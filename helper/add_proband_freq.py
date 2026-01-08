import pandas as pd
import numpy as np
import sqlite3 as sq
from sqlalchemy import create_engine


def add_freq(df):


    # Calculate total number of probands and non-probands
    total_proband_pt = df[df['is_proband'] == True]['pt_id'].nunique()
    total_nonproband_pt = df[df['is_proband'] == False]['pt_id'].nunique()

    # Cluster-wise unique pt_id count for probands only
    print('Getting proband-only cluster stats')
     # Step 1: Compute proband counts per cluster (excluding noise)
    proband_counts = (
        df[(df['is_proband'] == 1) & (df['cluster'] != -1)]
        .groupby('cluster')['pt_id']
        .nunique()
        .to_dict()
    )
    # Step 2: Assign row-level proband count correctly
    def get_proband_count(row):
        if row['cluster'] == -1:
            return 1 if row['is_proband'] else 0
        return proband_counts.get(row['cluster'], 0)

    df['proband_only_count'] = df.apply(get_proband_count, axis=1)
    df['proband_only_propensity'] = df['proband_only_count'] / total_proband_pt

    df['nonproband_only_count'] = df['unique_pt_id_count'] - df['proband_only_count']
    df['nonproband_only_propensity'] = df['nonproband_only_count'] / total_nonproband_pt

    df['UUID'] = df['chrom1'].astype(str)+'_'+df['pos1'].astype(str)+'_'+df['chrom2'].astype(str)+'_'+df['pos2'].astype(str)+'_'+df['SV_id'].astype(str)+'_'+df['pt_id'].astype(str)
    return df

def main(fn, table, db):
    df = pd.read_csv(f'./data/3.0.0/{fn}', sep='\t', index_col=False)
    print(df)
    df = add_freq(df)
    
    engine = create_engine(f'sqlite:///{db}.sqlite')  

    # Assuming you have 'row_id' in df and in the database
    update_cols = ['UUID', 'proband_only_count', 'proband_only_propensity', 'nonproband_only_count', 'nonproband_only_propensity']
    update_df = df[update_cols].copy()

    # Load existing table into a temporary table in the DB
    update_df.to_sql('tmp_update', engine, if_exists='replace', index=False)

    # Now run SQL to update original table using tmp_update
    conn = sq.connect(f'{db}.sqlite')
    cursor = conn.cursor()  
    cursor.execute(f"ALTER TABLE {table} ADD proband_only_count INT")
    cursor.execute(f"ALTER TABLE {table} ADD nonproband_only_count INT")
    cursor.execute(f"ALTER TABLE {table} ADD proband_only_propensity FLOAT")
    cursor.execute(f"ALTER TABLE {table} ADD nonproband_only_propensity FLOAT")
    cursor.execute(f"""
        UPDATE {table}
        SET
            proband_only_count = tmp.proband_only_count,
            nonproband_only_count = tmp.nonproband_only_count,
            proband_only_propensity = tmp.proband_only_propensity,
            nonproband_only_propensity = tmp.nonproband_only_propensity
        
        FROM tmp_update AS tmp
        WHERE {table}.UUID = tmp.UUID
    """)
    cursor.execute("""
        DROP TABLE tmp_update;
        """)

    conn.commit()
    conn.close()
   


if __name__ == '__main__':
    main('CLDB_P2_hg38_072425.tsv', 'P2_hg38', 'prod_CLDB_SR')