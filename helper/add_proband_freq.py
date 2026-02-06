import pandas as pd
import sqlite3 as sq
from sqlalchemy import create_engine
from dataclasses import dataclass

@dataclass
class UpdateJob:
    table_name: str
    input_file: str
    db_name: str

def calculate_frequencies(df: pd.DataFrame) -> pd.DataFrame:
    """Vectorized calculation of proband/non-proband stats."""
    print("Calculating frequencies...")
    
    # 1. Calculate Totals
    total_proband_pt = df.loc[df['IS_PROBAND'] == True, 'PT_ID'].nunique()
    total_nonproband_pt = df.loc[df['IS_PROBAND'] == False, 'PT_ID'].nunique()

    # 2. Vectorized Proband Counts per Cluster
    # Filter out noise (CLUSTER_ID -1) and non-probands
    mask = (df['IS_PROBAND'] == 1) & (df['CLUSTER_ID'] != -1)
    proband_counts = df[mask].groupby('CLUSTER_ID')['PT_ID'].nunique()

    # 3. Map counts back to DF (Handling Cluster -1 separately)
    df['proband_only_count'] = df['CLUSTER_ID'].map(proband_counts)
    
    # Fill noise/singular cases: If -1, use IS_PROBAND as 1/0
    noise_mask = df['CLUSTER_ID'] == -1
    df.loc[noise_mask, 'proband_only_count'] = df['IS_PROBAND'].astype(int)
    # Fill any clusters that had 0 probands
    df['proband_only_count'] = df['proband_only_count'].fillna(0).astype(int)

    # 4. Final Props
    df['proband_only_propensity'] = df['proband_only_count'] / total_proband_pt
    df['nonproband_only_count'] = df['UNIQUE_PT_COUNT'] - df['proband_only_count']
    df['nonproband_only_propensity'] = df['nonproband_only_count'] / total_nonproband_pt

    # 5. UUID Creation
    cols = ['chrom1', 'pos1', 'chrom2', 'pos2', 'SV_ID', 'PT_ID']
    df['UUID'] = df[cols].astype(str).agg('_'.join, axis=1)
    
    return df

def run_update_job(job: UpdateJob):
    print(f"--- Updating {job.table_name} from {job.input_file} ---")
    
    # Load and Process
    df = pd.read_csv(job.input_file, sep='\t', low_memory=False)
    df = calculate_frequencies(df)
    
    update_cols = ['UUID', 'proband_only_count', 'proband_only_propensity', 
                   'nonproband_only_count', 'nonproband_only_propensity']
    update_df = df[update_cols]

    # Database Operations
    db_path = job.db_name if job.db_name.endswith('.sqlite') else f"{job.db_name}.sqlite"
    engine = create_engine(f'sqlite:///{db_path}')
    
    with sq.connect(db_path) as conn:
        cursor = conn.cursor()
        
        # 1. Add columns (Ignore error if they already exist)
        new_cols = {
            "proband_only_count": "INT",
            "nonproband_only_count": "INT",
            "proband_only_propensity": "FLOAT",
            "nonproband_only_propensity": "FLOAT"
        }
        
        for col, dtype in new_cols.items():
            try:
                cursor.execute(f"ALTER TABLE {job.table_name} ADD COLUMN {col} {dtype}")
            except sq.OperationalError:
                pass # Column already exists
        
        # 2. Use Temporary Table for high-speed update
        update_df.to_sql('tmp_update', engine, if_exists='replace', index=False)
        
        # 3. Perform the Update
        # Note: Added an INDEX on tmp_update(UUID) to make the update instantaneous
        cursor.execute("CREATE INDEX idx_tmp_uuid ON tmp_update(UUID)")
        
        sql_update = f"""
            UPDATE {job.table_name}
            SET
                proband_only_count = (SELECT proband_only_count FROM tmp_update WHERE tmp_update.UUID = {job.table_name}.UUID),
                nonproband_only_count = (SELECT nonproband_only_count FROM tmp_update WHERE tmp_update.UUID = {job.table_name}.UUID),
                proband_only_propensity = (SELECT proband_only_propensity FROM tmp_update WHERE tmp_update.UUID = {job.table_name}.UUID),
                nonproband_only_propensity = (SELECT nonproband_only_propensity FROM tmp_update WHERE tmp_update.UUID = {job.table_name}.UUID)
            WHERE EXISTS (
                SELECT 1 FROM tmp_update WHERE tmp_update.UUID = {job.table_name}.UUID
            )
        """
        cursor.execute(sql_update)
        cursor.execute("DROP TABLE tmp_update")
        conn.commit()
    
    print(f"Update complete for {job.table_name}.")

# Example usage with multiple jobs
if __name__ == '__main__':
    jobs = [
        UpdateJob('sample_table', '/path/to/annotated/file.tsv', 'dev_sample_DB'),
        
        # Add more jobs here
    ]
    
    for job in jobs:
        run_update_job(job)
