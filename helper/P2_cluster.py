import gzip
import os
import pandas as pd
import numpy as np
import polars as pl
import argparse
from sklearn.cluster import DBSCAN


pd.set_option('display.expand_frame_repr', False)


def get_P2_df(meta_file, ref_ver):
	#global value for DBSCAN
	if ref_ver == 'hg19':
		## hg19
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
		# #hg38
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


	row_list = []
	meta = pd.read_csv(meta_file, sep='\t', index_col=False)
	print(f'Importing {meta.P2_path.nunique()} subjects')
	for m in meta.itertuples():
		pt_id = m.pt_id
		family = m.family
		project = m.project
		is_proband = m.is_proband
		system=m.system
		fname = m.P2_path
		if pd.isna(fname):
			continue
		print(f'{m.pt_id} from {m.project}')
		with gzip.open(fname, 'rt') as f:
			for line in f:
				tmp_dict={}
				if line.startswith('#'):
					continue
				else: 
					tmp = line.strip().split('\t')
					chrom1 = tmp[0][3:]
					pos1 = int(tmp[1])
					SV_id = tmp[2]
					tmp2 = tmp[7].split(';')
					SV_type = tmp2[3][7:]
					SV_len = round(float(tmp2[2][7:]))
					chrom2 = tmp2[5][5:].split('chr')[1] 
					pos2 = int(tmp2[6].split('=')[-1])
					genotype = tmp[-1].split(':')[0]

					#skipping for P2 artifact in chr2
					if ref_ver == 'hg19':
						#hg19
						if (chrom1 == '2') and (chrom2 == '2') and (33141211 <= pos1 <= 33141696):
							continue
						if (chrom1 == '2') and (chrom2 == '2') and (33141211 <= pos2 <= 33141696):
							continue
					elif ref_ver == 'hg38':
						# hg38 chr2:32916144-32916629
						if (chrom1 == '2') and (chrom2 == '2') and (32916144 <= pos1 <= 32916629):
							continue
						if (chrom1 == '2') and (chrom2 == '2') and (32916144 <= pos2 <= 32916629):
							continue

					tmp_dict.update({'chrom1':chrom1, 'pos1': pos1, 'chrom2': chrom2, 'pos2': pos2, 'SV_id': SV_id, 'SV_type': SV_type, 'SV_len': SV_len, 'genotype': genotype, 'pt_id': pt_id, 'family': family, 'project': project, 'is_proband': is_proband, 'system': system})
					row_list.append(tmp_dict)
	#load P2 files into single table
	df = pd.DataFrame(row_list)
	print(f'Imported {meta.P2_path.nunique()} subjects')



	#for testing
	# df = df[(df['chrom1'] == '1') & (df['pos1'] < 200000000)]
	# print(df.shape)

	#get cum position for DBSCAN
	df=df.merge(size_df, left_on='chrom1', right_on='Chromosome')
	df=df.drop(['Chromosome', 'Size'], axis = 1)
	df.rename(columns={df.columns[-1]: 'cum_size1'}, inplace=True)
	df=df.merge(size_df, left_on='chrom2', right_on='Chromosome')
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
		dbscan = DBSCAN(eps=epsilon, min_samples=min_samples, n_jobs=-2)
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



def main():
	P2 = get_P2_df('Z:/Members/clun/CLDB/meta/meta_SR_hg38.tsv', 'hg38')
	print('Writing Results to TSV')
	P2.to_csv('CLDB_P2_hg38_072425.tsv', sep='\t', index=False, chunksize=250000)
if __name__ == '__main__':
	main()