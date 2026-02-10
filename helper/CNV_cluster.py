import pandas as pd
import numpy as np
import os
import glob
from util.slmseg import * 
import rle
from joblib import Parallel, delayed
from tqdm import tqdm
from sklearn.cluster import DBSCAN

# Constants for Stage 2
REF_VER = 'hg19'
INPUT_META = 'Z:/Members/clun/CLDB/meta/meta_SR_hg19.tsv'
OUTPUT_DIR = 'Z:/Members/clun/CLDB/data/4.0.0/CLDB_CNV_hg19_intermediate/'
FINAL_FILE = 'Z:/Members/clun/CLDB/CLDB_CNV_hg19_021026.tsv'



def remove_chr(value):
    if value.startswith("chr"):
        return value.split('chr')[1]
    else:
        return value

def modify_chr(value):
	if value.startswith("chr"):
		return value
	else:
		return "chr" + value

def prep_df(df):
	df.columns = ['chr', 'start', 'end', 'coverage']
	df.iloc[:, 0] = df.iloc[:, 0].astype(str)
	df['chr'] = df['chr'].apply(modify_chr)
	return df

def normalization(rd):
	'''
	Normalize by chr
	'''
	rd['ratio'] = rd.groupby('chr')['coverage'].transform(lambda x: x / (x.median() + 0.00001))
	return rd

def get_seg(df, idx):
	'''
	get segments for single chromosome
	'''
	df = df[df['chr'] == idx]
	pos_ = df['start']
	tmp = df['ratio'] + 0.00001
	signal_ = np.log2(tmp)
	reset_seg_param()
	load_data(signal_, pos_)
	set_variables(0.3, 0.00001, 1000000, 0)
	SLM()
	results = data_seg()

	rle_ = rle.encode(results)
	end_indices = [sum(rle_[1][:i]) for i in range(1, len(rle_[1]) + 1)]
	out = pd.DataFrame({
		'chr': idx,
		'start': [0]+end_indices[:-1],
		'end': end_indices,
		'log2r': rle_[0]
		})
	out['start'] *=1000
	out['end'] *= 1000
	out['len'] = out['end'] - out['start']

	return out

def get_all_seg(df):
	'''
	get seg for all chromosomes
	'''
	out = []
	chromosomes = [f"chr{i}" for i in range(1, 23)] + ['chrX']  ##change this before deploy
	for c in chromosomes:
		tmp = get_seg(df, c)
		out.append(tmp)
	return pd.concat(out, ignore_index=True)

def mutate_log_lvl(log2r):
    conditions = [
        (log2r <= -1.525),
        ((log2r >= np.log2(1*0.9/2)) & (log2r <= np.log2(1*1.1/2))),
        ((log2r >= np.log2(2*0.9/2)) & (log2r <= np.log2(2*1.1/2))),
        ((log2r >= np.log2(3*0.9/2)) & (log2r <= np.log2(3*1.1/2))),
        ((log2r >= np.log2(4*0.9/2)) & (log2r <= np.log2(4*1.1/2))),
        (log2r >= 1.175)
    ]

    choices = ["HOM_DEL", "HET_DEL", "NML", "DUP", "TRP", "MUL_GAIN"]

    return np.select(conditions, choices, default="UND")

def get_all_cnv(df):
	'''
	wrapper to get all gain/loss cnv calls
	'''
	df = prep_df(df)
	x = normalization(df)
	y = get_all_seg(x)

	tmp=[]
	dup_ = y[(y['log2r'] > 0.4) & (y['log2r'] < 3) & (y['len'] >= 10000)].copy()
	dup_['type'] = 'gain'
	del_ = y[(y['log2r'] < -0.8) & (y['log2r'] > -5) & (y['len'] >= 10000)].copy()
	del_['type'] = 'loss'
	tmp.append(dup_)
	tmp.append(del_)
	out=pd.concat(tmp, ignore_index=True)
	out["SV_id"] = out.groupby("type").cumcount() + 1
	out["SV_id"] = out["SV_id"].astype(str).str.zfill(5)
	out["SV_id"] = out["type"] + "_" + out["SV_id"]
	out['cnv_lvl'] = mutate_log_lvl(out['log2r'])
	return out
# --- Utility Functions (Kept original logic, improved efficiency) ---

def mutate_log_lvl(log2r):
    conditions = [
        (log2r <= -1.525),
        ((log2r >= np.log2(0.9/2)) & (log2r <= np.log2(1.1/2))),
        ((log2r >= np.log2(1.8/2)) & (log2r <= np.log2(2.2/2))),
        ((log2r >= np.log2(2.7/2)) & (log2r <= np.log2(3.3/2))),
        ((log2r >= np.log2(3.6/2)) & (log2r <= np.log2(4.4/2))),
        (log2r >= 1.175)
    ]
    choices = ["HOM_DEL", "HET_DEL", "NML", "DUP", "TRP", "MUL_GAIN"]
    return np.select(conditions, choices, default="UND")

# --- Parallel Worker Task ---

def process_single_patient(m, output_dir):
    """Worker function: processes one patient and saves to disk."""
    try:
        # Check if already processed
        out_path = os.path.join(output_dir, f"{m.pt_id}.tsv")
        if os.path.exists(out_path):
            return True

        # Load and Prep
        df = pd.read_csv(m.MD_path, compression='gzip', sep='\t', dtype={0: str})
        
        # Core logic from your original script
        out = get_all_cnv(df) # Logic from your get_all_cnv
        
        # Post-processing
        out['pt_id'] = m.pt_id
        out['family'] = m.family
        out['project'] = m.project
        out['is_proband'] = m.is_proband
        out['system'] = m.system
        
        # Rename and Column Selection
        out['chr'] = out['chr'].apply(lambda x: x.split('chr')[1] if x.startswith("chr") else x)
        out['chrom1'] = out['chr']
        out['chrom2'] = out['chr']
        out = out.rename(columns={"start": "pos1", "end": "pos2", "cnv_lvl": "SV_type"})
        
        final_cols = ['chrom1', 'pos1', 'chrom2', 'pos2', 'SV_id', 'SV_type', 
                      'len', 'log2r', 'pt_id', 'family', 'project', 'is_proband', 'system']
        out = out[final_cols]
        
        
        out.to_csv(out_path, sep='\t', index=False)
        return True
    except Exception as e:
        print(f"Error on {m.pt_id}: {e}")
        return False

# --- Global Clustering Function ---

def apply_DBSCAN_refactored(df, ref_ver):
    """Optimized DBSCAN with vectorized cumulative position calculation."""
    chrom_configs={
		'hg19': {	'1':249250621 ,	'2':243199373 ,	'3':198022430 ,	'4':191154276 ,
					'5':180915260,	'6':171115067 ,	'7':159138663 ,	'8':146364022 ,
					'9':141213431 ,	'10':135534747,	'11':135006516,	'12':133851895,
					'13':115169878,	'14':107349540,	'15':102531392,	'16':90354753 ,
					'17':81195210 ,	'18':78077248 ,	'20':63025520 ,	'19':59128983 ,
					'22':51304566 ,	'21':48129895,	'X':155270560,	'Y':59373566},
		'hg38':{	'1'	:248956422,	'2'	:242193529,	'3'	:198295559,	'4'	:190214555,
					'5'	:181538259,	'6'	:170805979,	'7'	:159345973,	'X'	:156040895,	
					'8'	:145138636,	'9'	:138394717,	'11':135086622,	'10':133797422,	
					'12':133275309,	'13':114364328,	'14':107043718,	'15':101991189,	
					'16':90338345,	'17':83257441,	'18':80373285,	'20':64444167,	
					'19':58617616,	'Y'	:57227415,	'22':50818468,	'21':46709983}
	}
    # Choose sizes
    sizes = chrom_configs[ref_ver] 
    
    # Vectorized Cumulative Position
    size_df = pd.DataFrame(list(sizes.items()), columns=['chrom1', 'Size'])
    size_df['cum_offset'] = size_df['Size'].cumsum().shift(fill_value=0)
    
    df = df.merge(size_df[['chrom1', 'cum_offset']], on='chrom1', how='left')
    df['cum_pos1'] = df['pos1'] + df['cum_offset']
    df['cum_pos2'] = df['pos2'] + df['cum_offset']

    # DBSCAN
    df['cluster'] = -1
    cluster_counter = 0
    
    for sv_type, group in df.groupby('SV_type'):
        if len(group) < 2: continue
        
        coords = group[['cum_pos1', 'cum_pos2']]
        db = DBSCAN(eps=500, min_samples=2).fit(coords)
        
        # Global unique cluster IDs
        labels = db.labels_
        mask = labels != -1
        labels[mask] = labels[mask] + cluster_counter + 1
        df.loc[group.index, 'cluster'] = labels
        
        if any(mask):
            cluster_counter = labels.max()

    # Calculate frequencies (Vectorized)
    n_sub = df['pt_id'].nunique()
    df['count'] = df.groupby('cluster')['cluster'].transform('count')
    df.loc[df['cluster'] == -1, 'count'] = 1

    df['unique_pt_id_count'] = df.groupby('cluster')['pt_id'].transform('nunique')
    df.loc[df['cluster'] == -1, 'unique_pt_id_count'] = 1
    df['pseudo_db_freq'] = df['unique_pt_id_count'] / n_sub

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
    df=df.drop(['cum_offset', 'cum_pos2', 'cum_pos1'], axis = 1)

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

# --- Main Execution ---

def main(meta_file, output_dir, final_output, ref_ver):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    meta = pd.read_csv(meta_file, sep='\t')
    meta = meta[meta["CNV_outlier"] == 0]
    
    # STAGE 1: Parallel Process Patients
    print(f"Stage 1: Processing {len(meta)} patients in parallel...")
    Parallel(n_jobs=-1)(
        delayed(process_single_patient)(m, output_dir) 
        for m in tqdm(list(meta.itertuples()), desc="Processing CNVs")
    )

    # STAGE 2: Bulk Load and Cluster
    print("Stage 2: Bulk loading intermediate files...")
    all_files = glob.glob(os.path.join(output_dir, "*.tsv"))
    # Fast load: using list comprehension + concat
    df_all = pd.concat([pd.read_csv(f, sep='\t', dtype={'chrom1': str}) for f in all_files], ignore_index=True)
    
    print("Applying DBSCAN Clustering...")
    df_clustered = apply_DBSCAN_refactored(df_all, ref_ver)

    print(f"Saving final result to {final_output}")
    df_clustered.to_csv(final_output, sep='\t', index=False)

if __name__ == '__main__':
    main(INPUT_META, 
         OUTPUT_DIR, 
         FINAL_FILE, 
         REF_VER)