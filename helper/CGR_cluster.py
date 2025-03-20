import pandas as pd
import numpy as np
import sqlite3 as sq
from sqlalchemy import create_engine

pd.set_option('display.expand_frame_repr', False)


def modify_chr(value):
    if value.startswith("chr"):
        return value
    else:
        return "chr" + value

        
def query(pt_id):
    conn = sq.connect(f'GREGoR_052224.sqlite')
    cnv = pd.read_sql(f'select chr, start, end, len, type, cnv_lvl, pt_id, family, is_proband, project, cluster, count, unique_pt_id_count, psuedo_df_freq, SD_Overlap, OMIM_Count, RefSeq_Count from CNV_hg38 where 1=1 AND pt_id = "{pt_id}"', conn)
    p2 = pd.read_sql(f'select chrom1 AS chr, pos1 AS start, pos2 AS end, SV_id, SV_type from P2_hg38 where 1=1 AND pt_id = "{pt_id}" AND SV_type != "BND"', conn)
    p2['chr']=p2['chr'].astype(str)
    p2['chr']=p2['chr'].apply(modify_chr)
    conn.close()
    return p2, cnv


def process_cnv(p2, cnv):
    # Filter p2
    p2 = p2[abs(p2['start'] - p2['end']) > 1000]
    p2 = p2.reset_index(drop=True)
    # print(p2)
    p2_dict={}
    for name, group in p2.groupby('chr'):
        p2_dict[name] = group
    # Process start positions
    final_idx=[]
    for idx, row in cnv.iterrows():
        _chr=row['chr']
        abs_diff_start = abs(row['start'] - p2_dict[_chr]['start'])
        closest_idx_start = abs_diff_start.idxmin()
        abs_diff_end = abs(row['start'] - p2_dict[_chr]['end'])
        closest_idx_end = abs_diff_end.idxmin()
        if  abs_diff_start[closest_idx_start] < abs_diff_end[closest_idx_end]:
            final_idx.append(closest_idx_start)
        else:
            final_idx.append(closest_idx_end)
        
    left_df=p2.iloc[final_idx,]
    left_df.columns='left_'+left_df.columns
    left_df['left_width'] = abs(left_df['left_start']- left_df['left_end'])
    left_df=left_df.reset_index(drop=True)
    
    # Process end positions    
    final_idx=[]
    for idx, row in cnv.iterrows():
        _chr=row['chr']
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
    right_df['right_width'] = abs(right_df['right_start']- right_df['right_end'])    
    right_df=right_df.reset_index(drop=True)

    merged_df = pd.concat([cnv,left_df, right_df], axis=1)
    
    # Calculate distances
    merged_df['left_dist'] = merged_df.apply(lambda x: min(abs(x['start']-x['left_start']), abs(x['start']-x['left_end'])), axis=1)
    merged_df['right_dist'] = merged_df.apply(lambda x: min(abs(x['start']-x['right_start']), abs(x['start']-x['right_end'])), axis=1)

    # Matching criteria
    merged_df['match_SV'] = (merged_df['left_SV_id'] == merged_df['right_SV_id']).astype(int)
    merged_df['match_left_SV_size'] = np.where(abs(merged_df['len'] - merged_df['left_width']) < 2000, 1,0)
    merged_df['match_right_SV_size'] = (abs(merged_df['len'] - merged_df['right_width']) < 2000).astype(int)
    merged_df['match_left_SV_type'] = ((merged_df['type'] == 'gain') & (merged_df['left_SV_type'] == 'DUP') | (merged_df['type'] == 'loss') & (merged_df['left_SV_type'] == 'DEL')).astype(int)
    merged_df['match_right_SV_type'] = ((merged_df['type'] == 'gain') & (merged_df['right_SV_type'] == 'DUP') | (merged_df['type'] == 'loss') & (merged_df['right_SV_type'] == 'DEL')).astype(int)
    merged_df['match_CNV_SV'] = ((merged_df['match_SV'] == 1) & (merged_df['match_left_SV_type'] == 1) & (merged_df['match_left_SV_size'] == 1)).astype(int)
    merged_df = merged_df.sort_values(by='chr')
    return merged_df


def main():
    df = pd.read_csv('meta/gregor_meta.tsv', sep='\t')
    df = df[['pt_id', 'P2_path']]
    df=df.dropna().reset_index()
    pt_ids=df['pt_id']
    for pt_id in pt_ids:
        print(pt_id)
        p2, cnv = query(pt_id)
        out = process_cnv(p2, cnv)
        out.to_csv(f'CGR_gregor/{pt_id}.tsv', sep='\t', index=False)

if __name__ == '__main__':
    main()