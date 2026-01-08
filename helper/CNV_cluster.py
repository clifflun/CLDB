import pandas as pd
import numpy as np
from util.slmseg import * 
from itertools import groupby
import rle
from sklearn.cluster import DBSCAN


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

def postproc(df, **kwargs):
	df['pt_id'] = kwargs['pt_id']
	df['family'] = kwargs['family']
	df['is_proband'] = kwargs['is_proband']
	df['project'] = kwargs['project']
	df['system'] = kwargs['system']

def apply_DBSCAN(df, ref_ver):
	#global value for DBSCAN
	# #hg19
	if ref_ver == 'hg19':
		chrom_size = {
			'1':249250621 ,
			'2':243199373 ,
			'3':198022430 ,
			'4':191154276 ,
			'5':180915260,
			'6':171115067 ,
			'7':159138663 ,
			'8':146364022 ,
			'9':141213431 ,
			'10':135534747,
			'11':135006516,
			'12':133851895,
			'13':115169878,
			'14':107349540,
			'15':102531392,
			'16':90354753 ,
			'17':81195210 ,
			'18':78077248 ,
			'20':63025520 ,
			'19':59128983 ,
			'22':51304566 ,
			'21':48129895,
			'X':155270560,
			'Y':59373566 }
	elif ref_ver == 'hg38':
		#hg38
		chrom_size={
			'1'	:248956422,
			'2'	:242193529,
			'3'	:198295559,
			'4'	:190214555,
			'5'	:181538259,
			'6'	:170805979,
			'7'	:159345973,
			'X'	:156040895,
			'8'	:145138636,
			'9'	:138394717,
			'11'	:135086622,
			'10'	:133797422,
			'12'	:133275309,
			'13'	:114364328,
			'14'	:107043718,
			'15'	:101991189,
			'16'	:90338345,
			'17'	:83257441,
			'18'	:80373285,
			'20'	:64444167,
			'19'	:58617616,
			'Y'	:57227415,
			'22'	:50818468,
			'21'	:46709983
		}

	size_df = pd.DataFrame(list(chrom_size.items()), columns=['Chromosome', 'Size'])
	size_df['cum_size'] = size_df['Size'].cumsum().shift(fill_value=0)


	#get cum position for DBSCAN
	df=df.merge(size_df, left_on='chrom1', right_on='Chromosome')
	print(df)
	
	df=df.drop(['Chromosome', 'Size'], axis = 1)
	df.rename(columns={df.columns[-1]: 'cum_size1'}, inplace=True)
	df=df.merge(size_df, left_on='chrom1', right_on='Chromosome')
	df.rename(columns={df.columns[-1]: 'cum_size2'}, inplace=True)

	df['cum_pos1'] = df['pos1']+df['cum_size1']
	df['cum_pos2'] = df['pos2']+df['cum_size2']

	#DBSCAN parameters
	epsilon = 500
	min_samples = 2
	cluster_counter = 0

	#DBSCAN by SV type
	grouped_SV_type = df.groupby('SV_type')

	for name, group in grouped_SV_type:
		print(f'Clustering {name}')
		features = ['cum_pos1', 'cum_pos2']
		X = group[features]
		dbscan = DBSCAN(eps=epsilon, min_samples=min_samples)
		group['cluster_tmp'] = dbscan.fit_predict(X)
		group['cluster'] = np.where(group['cluster_tmp'] != -1, group['cluster_tmp']+cluster_counter+1, -1)
		cluster_counter = group['cluster'].max()
		df.loc[group.index, 'cluster'] = group['cluster']

	df=df.drop(['Chromosome', 'Size', 'cum_size1', 'cum_size2', 'cum_pos2', 'cum_pos1'], axis = 1)
	total_pt = len(df['pt_id'].unique())

	#get count and propensity
	print('Getting cluster count')
	df['count'] = df.groupby(['cluster'])['cluster'].transform('count')
	df.loc[df['cluster'] == -1, 'count'] = 1
	df['cluster_propensity'] = df['count']/total_pt
	
	#get pseudo db freq
	print('Getting pseudo db freq')
	df['unique_pt_id_count'] = df.groupby('cluster')['pt_id'].transform('nunique')
	df['unique_pt_id_count'] = np.where(df['cluster'] == -1, 1, df['unique_pt_id_count'])
	nsub=df.pt_id.nunique()
	df['psuedo_df_freq'] = df['unique_pt_id_count']/nsub
	return df

def collapse(df):
	print('Collapsing calls in repeats')
	grouped = df.groupby(['type', 'IDR_Disrupt_left', 'IDR_Disrupt_right'])['cluster'].agg(['unique'])
	grouped = grouped[grouped['unique'].apply(len) > 1]
	for row in grouped.itertuples():
		#avoid grouping calls without IDR disrupt
		if row[0][1] == '' and row[0][1] == '':
			continue
		print(row)
		idx_ = df[(df['type'] == row[0][0]) & (df['IDR_Disrupt_left'] == row[0][1]) & (df['IDR_Disrupt_right'] == row[0][2])].index
		new_id = df.loc[idx_]['cluster'].max()
		df.loc[idx_, 'cluster'] = new_id


	##recalculate freq
	#get count and propensity
	print('Getting cluster count')
	total_pt = len(df['pt_id'].unique())
	df['count'] = df.groupby(['cluster'])['cluster'].transform('count')
	df.loc[df['cluster'] == -1, 'count'] = 1
	df['cluster_propensity'] = df['count']/total_pt
	
	#get pseudo db freq
	print('Getting pseudo db freq')
	df['unique_pt_id_count'] = df.groupby('cluster')['pt_id'].transform('nunique')
	df['unique_pt_id_count'] = np.where(df['cluster'] == -1, 1, df['unique_pt_id_count'])
	nsub=df.pt_id.nunique()
	df['psuedo_df_freq'] = df['unique_pt_id_count']/nsub
	return df

def main(fn, fo, ref_ver):
	tmp=[]
	meta = pd.read_csv(fn, sep='\t', index_col=False)
	meta = meta[meta["CNV_outlier"] == 0]
	for m in meta.itertuples():
		fname=m.MD_path 
		pt_id=m.pt_id
		family=m.family
		project=m.project
		is_proband=m.is_proband
		system=m.system
		print(f'{pt_id} from {project}')
		df = pd.read_csv(fname, compression='gzip', sep='\t', dtype={0: str})
		out = get_all_cnv(df)
		postproc(out, pt_id=pt_id, family=family, project=project, is_proband=is_proband, system=system)
		out['chr'] = out['chr'].apply(remove_chr)
		out['chrom2'] = out['chr']
		out = out.rename(columns={"start": "pos1", "end": "pos2", "chr": "chrom1", "cnv_lvl": "SV_type"})
		out = out[['chrom1', 'pos1', 'chrom2', 'pos2', 'SV_id', 'SV_type', 'len', 'log2r', 'pt_id', 'family', 'project', 'is_proband', 'system']]
		tmp.append(out)
	df = pd.concat(tmp, ignore_index=True)
	# df=pd.read_csv(f'{fo}.tsv', sep='\t', index_col=False, dtype ={'chrom1': str})
	df = apply_DBSCAN(df, ref_ver)
	df.to_csv(f'{fo}.tsv', sep='\t', index=False)

	

if __name__ == '__main__':
	pd.set_option('display.expand_frame_repr', False)
	main('Z:/Members/clun/CLDB/meta/meta_SR_hg38.tsv', 'Z:/Members/clun/CLDB/CLDB_CNV_hg38_062625', 'hg38')