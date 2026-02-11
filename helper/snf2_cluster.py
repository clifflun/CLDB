import gzip
import os
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
import re
from joblib import Parallel, delayed
from tqdm import tqdm

pd.set_option('display.expand_frame_repr', False)

# --- CONFIGURATION ---
REF_VER = 'hg38' 
META_FILE = 'Z:/Members/clun/CLDB/meta/meta_LR.tsv'
OUTPUT_FILE = f'CLDB_snf2_{REF_VER}.tsv'

# ---------------------------------------------------------------------------
# Worker Task: Parallel Parsing of gzipped files
# ---------------------------------------------------------------------------

def parse_single_vcf(m, ref_ver):
    """
    Worker function to parse a single subject's file.
    """
    row_list = []
    pt_id, family, project, is_proband = m.pt_id, m.family, m.project, m.is_proband
    
    if ref_ver == 'hg38':
        fname, system = m.hg38_path, 'Revio'
    elif ref_ver == 'hg19':
        fname, system = m.hg19_path, 'SII'
    else:
        return []

    if pd.isna(fname) or not os.path.exists(fname):
        return []

    try:
        with gzip.open(fname, 'rt') as f:
            for line in f:
                if line.startswith('#'): continue
                
                parts = line.strip().split('\t')
                chrom1 = parts[0].replace('chr', '')
                pos1 = int(parts[1])
                SV_id = parts[2]
                info_raw = parts[7]
                
                # Faster info parsing
                info_dict = dict(item.split("=", 1) if "=" in item else (item, True) 
                                for item in info_raw.split(";"))
                
                SV_type = info_dict.get('SVTYPE', 'UNK')
                
                if SV_type != 'BND':
                    SV_len = abs(round(float(info_dict.get('SVLEN', 0))))
                    chrom2 = chrom1
                    pos2 = int(info_dict.get('END', pos1 + SV_len))
                else:
                    SV_len = 0
                    match = re.search(r'[\[\]]([^\[\]]+)[\[\]]', parts[4])
                    if match:
                        dest = match.group(1).split(':')
                        chrom2 = dest[0].replace('chr', '')
                        pos2 = int(dest[1])
                    else:
                        chrom2, pos2 = chrom1, pos1
                
                genotype = parts[-1].split(':')[0]
                
                row_list.append({
                    'chrom1': chrom1, 'pos1': pos1, 'chrom2': chrom2, 'pos2': pos2, 
                    'SV_id': SV_id, 'SV_type': SV_type, 'SV_len': SV_len, 
                    'genotype': genotype, 'pt_id': pt_id, 'family': family, 
                    'project': project, 'is_proband': is_proband, 'system': system
                })
    except Exception:
        return []
    
    return row_list

# ---------------------------------------------------------------------------
# Core Logic Functions
# ---------------------------------------------------------------------------

def get_chrom_sizes(ref_ver):
    # (Kept original logic for size dictionary)
    sizes = {
        'hg38': {'1': 248956422, '2': 242193529, '3': 198295559, '4': 190214555, '5': 181538259, '6': 170805979, '7': 159345973, 'X': 156040895, '8': 145138636, '9': 138394717, '11': 135086622, '10': 133797422, '12': 133275309, '13': 114364328, '14': 107043718, '15': 101991189, '16': 90338345, '17': 83257441, '18': 80373285, '20': 64444167, '19': 58617616, 'Y': 57227415, '22': 50818468, '21': 46709983},
        'hg19': {'1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276, '5': 180915260, '6': 171115067, '7': 159138663, '8': 146364022, '9': 141213431, '10': 135534747, '11': 135006516, '12': 133851895, '13': 115169878, '14': 107349540, '15': 102531392, '16': 90354753, '17': 81195210, '18': 78077248, '20': 63025520, '19': 59128983, '22': 51304566, '21': 48129895, 'X': 155270560, 'Y': 59373566}
    }[ref_ver]
    size_df = pd.DataFrame(list(sizes.items()), columns=['Chromosome', 'Size'])
    size_df['cum_size'] = size_df['Size'].cumsum().shift(fill_value=0)
    return size_df

def cluster(df, size_df):
    print('Clustering using DBSCAN')
    
    df = df.merge(size_df, left_on='chrom1', right_on='Chromosome')
    df = df.drop(['Chromosome', 'Size'], axis=1)
    df.rename(columns={df.columns[-1]: 'cum_size1'}, inplace=True)
    
    df = df.merge(size_df, left_on='chrom2', right_on='Chromosome')
    df.rename(columns={df.columns[-1]: 'cum_size2'}, inplace=True)

    df['cum_pos1'] = df['pos1'] + df['cum_size1']
    df['cum_pos2'] = df['pos2'] + df['cum_size2']

    epsilon = 500
    min_samples = 2
    cluster_counter = 0

    grouped_SV_type = df.groupby('SV_type')
    
    for name, group in grouped_SV_type:
        print(f'Clustering {name}')
        features = ['cum_pos1', 'cum_pos2']
        X = group[features]
        
        # Optimized: Parallel Processing
        dbscan = DBSCAN(eps=epsilon, min_samples=min_samples, n_jobs=-2)
        
        group['cluster_tmp'] = dbscan.fit_predict(X)
        group['cluster'] = np.where(group['cluster_tmp'] != -1, group['cluster_tmp'] + cluster_counter + 1, -1)
        
        if not group['cluster'].empty:
            cluster_counter = group['cluster'].max()
            
        df.loc[group.index, 'cluster'] = group['cluster']

    df = df.drop(['Chromosome', 'Size', 'cum_size1', 'cum_size2', 'cum_pos2', 'cum_pos1'], axis=1)

    total_pt = len(df['pt_id'].unique())
    print(f"SV Types processed: {df.SV_type.unique()}")

    print('Getting cluster count')
    df['count'] = df.groupby(['cluster'])['cluster'].transform('count')
    df.loc[df['cluster'] == -1, 'count'] = 1
    df['cluster_propensity'] = df['count'] / total_pt

    print('Getting pseudo db freq')
    df['unique_pt_id_count'] = df.groupby('cluster')['pt_id'].transform('nunique')
    df['unique_pt_id_count'] = np.where(df['cluster'] == -1, 1, df['unique_pt_id_count'])
    nsub = df.pt_id.nunique()
    df['psuedo_df_freq'] = df['unique_pt_id_count'] / nsub


    # 1. Calculate Totals
    total_proband_pt = df.loc[df['is_proband'] == True, 'pt_id'].nunique()
    total_nonproband_pt = df.loc[df['is_proband'] == False, 'pt_id'].nunique()

    # 2. Vectorized Proband Counts per Cluster
    # Filter out noise (CLUSTER_ID -1) and non-probands
    mask = (df['is_proband'] == 1) & (df['cluster'] != -1)
    proband_counts = df[mask].groupby('cluster')['pt_id'].nunique()
    df['proband_only_count'] = df['cluster'].map(proband_counts)

    # Fill noise/singular cases: If -1, use is_proband as 1/0
    noise_mask = df['cluster'] == -1
    df.loc[noise_mask, 'proband_only_count'] = df['is_proband'].astype(int)
    df['proband_only_count'] = df['proband_only_count'].fillna(0).astype(int)
    df['proband_only_propensity'] = df['proband_only_count'] / total_proband_pt
    df['nonproband_only_count'] = df['unique_pt_id_count'] - df['proband_only_count']
    df['nonproband_only_propensity'] = df['nonproband_only_count'] / total_nonproband_pt

    # --- FINAL RENAMING TO MATCH DB SCHEMA ---
    # Mixed Case Mapping based on verified DB schema
    rename_map = {
        'SV_id': 'SV_ID',
        'SV_type': 'SV_TYPE',
        'SV_len': 'SV_LEN',
        'pt_id': 'PT_ID',
        'family': 'FAM_ID',
        'project': 'PROJECT',
        'is_proband': 'IS_PROBAND',
        'system': 'SYSTEM',
        'cluster': 'CLUSTER_ID',
        'count': 'COUNT',
        'cluster_propensity': 'CLUSTER_PROPENSITY',
        'unique_pt_id_count': 'UNIQUE_PT_COUNT',
        'psuedo_df_freq': 'PSEUDO_FREQ'
    }
    df = df.rename(columns=rename_map)
    return df
def main():
    size_df = get_chrom_sizes(REF_VER)
    meta = pd.read_csv(META_FILE, sep='\t', low_memory=False)
    
    # --- PARALLEL EXECUTION ---
    print(f"Starting Parallel Parse for {len(meta)} subjects...")
    # n_jobs=-1 uses all cores. 
    # Use delayed() to pass each meta row to the parser.
    results = Parallel(n_jobs=-1)(
        delayed(parse_single_vcf)(m, REF_VER) 
        for m in tqdm(list(meta.itertuples()), desc="Parsing Files")
    )
    
    # Flatten the list of lists into a single DataFrame
    flat_results = [item for sublist in results for item in sublist]
    df = pd.DataFrame(flat_results)
    # df.to_csv('debug.tsv', sep='\t', index=None)
    if df.empty:
        print("No data found."); return

    df_clustered = cluster(df, size_df)
    
    print(f'Writing to {OUTPUT_FILE}...')
    # Dropping technical columns before save
    cols_to_drop = ['cum_pos1', 'cum_pos2']
    df_clustered.drop(columns=cols_to_drop, errors='ignore')

    print('converting to pl')
	pl_df = pl.from_pandas(df_clustered)
	print('writing output')
	pl_df.write_csv(OUTPUT_NAME, separator='\t')
    print('Done')

if __name__ == '__main__':
    main()