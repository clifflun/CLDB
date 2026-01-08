import pandas as pd
import numpy as np
import sqlite3 as sq
from sqlalchemy import create_engine
import os

pd.set_option('display.expand_frame_repr', False)

def modify_chr(value):
    if value.startswith("chr"):
        return value
    else:
        return "chr" + value

def query(pt_id, conn):
    # Use the passed connection instead of creating a new one every time
    cnv = pd.read_sql(f'SELECT UUID, chrom1 AS chr, pos1 AS start, chrom2, pos2 AS end, SV_ID AS CNV_ID, SV_TYPE AS CNV_TYPE, SV_LEN AS CNV_LEN, PT_ID, FAM_ID, IS_PROBAND, PROJECT, CLUSTER_ID, COUNT, UNIQUE_PT_COUNT, PSEUDO_FREQ, SD_overlap, OMIM_count, RefSeq_count from CNV_hg38 where 1=1 AND pt_id = "{pt_id}"', conn)
    p2 = pd.read_sql(f'SELECT chrom1 AS chr, pos1 AS start, pos2 AS end, SV_ID, SV_TYPE from P2_hg38 where 1=1 AND pt_id = "{pt_id}" AND SV_TYPE != "BND"', conn)
    p2['chr']=p2['chr'].astype(str)
    cnv['chr']=cnv['chr'].astype(str)
    return p2, cnv

def process_cnv(p2, cnv):
    # Filter p2
    p2 = p2[abs(p2['start'] - p2['end']) > 1000]
    p2 = p2.reset_index(drop=True)
    p2_dict={}
    for name, group in p2.groupby('chr'):
        p2_dict[name] = group
        
    # Process start positions
    final_idx=[]
    for idx, row in cnv.iterrows():
        _chr=row['chr']
        if _chr not in p2_dict:
            continue

        abs_diff_start = abs(row['start'] - p2_dict[_chr]['start'])
        closest_idx_start = abs_diff_start.idxmin()
        abs_diff_end = abs(row['start'] - p2_dict[_chr]['end'])
        closest_idx_end = abs_diff_end.idxmin()
        if  abs_diff_start[closest_idx_start] < abs_diff_end[closest_idx_end]:
            final_idx.append(closest_idx_start)
        else:
            final_idx.append(closest_idx_end)
        
    if not final_idx: 
        return pd.DataFrame()

    left_df=p2.iloc[final_idx,]
    left_df.columns='left_'+left_df.columns
    # Fix SettingWithCopyWarning by operating on the copy directly or creating a new column properly
    left_df = left_df.copy() # Ensure it's a copy, not a view
    left_df['left_width'] = abs(left_df['left_start']-left_df['left_end'])
    left_df=left_df.reset_index(drop=True)
    
    # Process end positions    
    final_idx=[]
    for idx, row in cnv.iterrows():
        _chr=row['chr']
        if _chr not in p2_dict:
            continue

        abs_diff_start = abs(row['end'] - p2_dict[_chr]['start'])
        closest_idx_start = abs_diff_start.idxmin()
        abs_diff_end = abs(row['end'] - p2_dict[_chr]['end'])
        closest_idx_end = abs_diff_end.idxmin()
        if  abs_diff_start[closest_idx_start] < abs_diff_end[closest_idx_end]:
            final_idx.append(closest_idx_start)
        else:
            final_idx.append(closest_idx_end)
        
    right_df=p2.iloc[final_idx,]
    right_df.columns='right_'+right_df.columns
    # Fix SettingWithCopyWarning
    right_df = right_df.copy()
    right_df['right_width'] = abs(right_df['right_start']- right_df['right_end'])    
    right_df=right_df.reset_index(drop=True)

    # Re-aligning CNV to match the filtered results if necessary
    merged_df = pd.concat([cnv.reset_index(drop=True), left_df, right_df], axis=1)
    
    # Calculate distances
    merged_df['left_dist'] = merged_df.apply(lambda x: min(abs(x['start']-x['left_start']), abs(x['start']-x['left_end'])), axis=1)
    merged_df['right_dist'] = merged_df.apply(lambda x: min(abs(x['start']-x['right_start']), abs(x['start']-x['right_end'])), axis=1)

    # Matching criteria
    merged_df['match_SV'] = (merged_df['left_SV_ID'] == merged_df['right_SV_ID']).astype(int)
    merged_df['match_left_SV_size'] = np.where(abs(merged_df['CNV_LEN'] - merged_df['left_width']) < 2000, 1,0)
    merged_df['match_right_SV_size'] = (abs(merged_df['CNV_LEN'] - merged_df['right_width']) < 2000).astype(int)
    merged_df['match_left_SV_TYPE'] = (((merged_df['CNV_TYPE'] == 'DUP') | (merged_df['CNV_TYPE'] == 'TRP')) & (merged_df['left_SV_TYPE'] == 'DUP')) | (((merged_df['CNV_TYPE'] == 'HET_DEL') | (merged_df['CNV_TYPE'] == 'HOM_DEL')) & (merged_df['left_SV_TYPE'] == 'DEL')).astype(int)
    merged_df['match_right_SV_TYPE'] = (((merged_df['CNV_TYPE'] == 'DUP') | (merged_df['CNV_TYPE'] == 'TRP')) & (merged_df['right_SV_TYPE'] == 'DUP')) | (((merged_df['CNV_TYPE'] == 'HET_DEL') | (merged_df['CNV_TYPE'] == 'HOM_DEL')) & (merged_df['right_SV_TYPE'] == 'DEL')).astype(int)
    merged_df['match_CNV_SV'] = ((merged_df['match_SV'] == 1) & (merged_df['match_left_SV_TYPE'] == 1) & (merged_df['match_left_SV_size'] == 1)).astype(int)
    merged_df = merged_df.sort_values(by='chr')
    merged_df.rename(columns={'chr':'chrom1', 'start':'pos1', 'end':'pos2', 'CNV_ID': 'SV_ID', 'CNV_TYPE': 'SV_TYPE', 'CNV_LEN': 'SV_LEN'}, inplace=True)
    return merged_df


def main():
    # Load metadata
    df = pd.read_csv('Z:/Members/clun/CLDB/meta/meta_SR_hg38.tsv', sep='\t')
    df = df[['pt_id', 'P2_path']]
    df=df.dropna().reset_index()
    pt_ids=df['pt_id']
    
    output_dir = './data/3.0.0/CLDB_CGR_hg38_123125'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")

    # OPEN DB CONNECTION ONCE
    db_path = 'prod_CLDB_SR.sqlite'
    conn = sq.connect(db_path)
    print(f"Connected to database: {db_path}")

    try:
        for pt_id in pt_ids:
            print(pt_id)
            try:
                # Pass connection to query function
                p2, cnv = query(pt_id, conn)
                
                if p2.empty or cnv.empty:
                    # print(f"  Skipping {pt_id}: No P2 or CNV data found.")
                    continue
                    
                out = process_cnv(p2, cnv)
                
                if out.empty:
                    continue
                    
                out.to_csv(f'{output_dir}/{pt_id}.tsv', sep='\t', index=False)
            except Exception as e:
                print(f"Error processing {pt_id}: {e}")
    finally:
        # Always close connection at the end
        conn.close()
        print("Database connection closed.")

if __name__ == '__main__':
    main()