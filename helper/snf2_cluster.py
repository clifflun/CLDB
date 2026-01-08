import gzip
import os
import pandas as pd
import numpy as np
import polars as pl
import argparse
import re
import glob
from sklearn.cluster import DBSCAN
import sqlite3 as sq
from sqlalchemy import create_engine

pd.set_option('display.expand_frame_repr', False)

#global value for DBSCAN
## hg19
# chrom_size = {
# 	'1':249250621 ,
# 	'2':243199373 ,
# 	'3':198022430 ,
# 	'4':191154276 ,
# 	'5':180915260,
# 	'6':171115067 ,
# 	'7':159138663 ,
# 	'8':146364022 ,
# 	'9':141213431 ,
# 	'10':135534747,
# 	'11':135006516,
# 	'12':133851895,
# 	'13':115169878,
# 	'14':107349540,
# 	'15':102531392,
# 	'16':90354753 ,
# 	'17':81195210 ,
# 	'18':78077248 ,
# 	'20':63025520 ,
# 	'19':59128983 ,
# 	'22':51304566 ,
# 	'21':48129895,
# 	'X':155270560,
# 	'Y':59373566 }

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

def load_df(meta_file):
	row_list = []
	meta = pd.read_csv(meta_file, sep='\t', index_col=False)
	print(f'Searching {meta.pt_id.nunique()} subjects')
	for m in meta.itertuples():
		pt_id = m.pt_id
		family = m.family
		project = m.project
		is_proband = m.is_proband

		# fname = m.hg19_path
		# system='SII'
		fname = m.hg38_path
		system='Revio'
		if pd.isna(fname):
			continue

		print(f'{m.pt_id} from {m.project}')
		with gzip.open(fname, 'rt') as f:
			for line in f:
				tmp_dict={}
				if line.startswith('#'):
					continue
				else: 
					line = line.strip().split('\t')
					chrom1 = line[0].replace('chr', '')
					pos1 = int(line[1])
					SV_id = line[2]
					info = line[7].split(';')

					SV_type = info[1][7:]
					if SV_type != 'BND':
						SV_len = abs(round(float(info[2][6:])))
						chrom2=chrom1
						pos2 = int(info[3].split('=')[-1])
					else:
						SV_len = 0
						chrom2 = info[6][8:].replace('chr', '') #change to 5: if no chr
						pos2 = int(re.split("[\[\]]", line[4].split(':')[1])[0])			
					genotype = line[-1].split(':')[0]
					
					tmp_dict.update({'chrom1':chrom1, 'pos1': pos1, 'chrom2': chrom2, 'pos2': pos2, 'SV_id': SV_id, 'SV_type': SV_type, 'SV_len': SV_len, 'genotype': genotype, 'pt_id': pt_id, 'family': family, 'project': project, 'is_proband': is_proband, 'system': system})
					row_list.append(tmp_dict)
	#load files into single table
	df = pd.DataFrame(row_list)
	print(f'Imported {df["pt_id"].nunique()} subjects')
	return df

def cluster(df):
	print('Clustering using DBSCAN')
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
		# print(group)
		features = ['cum_pos1', 'cum_pos2']
		X = group[features]
		dbscan = DBSCAN(eps=epsilon, min_samples=min_samples)
		group['cluster_tmp'] = dbscan.fit_predict(X)
		group['cluster'] = np.where(group['cluster_tmp'] != -1, group['cluster_tmp']+cluster_counter+1, -1)
		cluster_counter = group['cluster'].max()
		df.loc[group.index, 'cluster'] = group['cluster']

	df=df.drop(['Chromosome', 'Size', 'cum_size1', 'cum_size2', 'cum_pos2', 'cum_pos1'], axis = 1)
	total_pt = len(df['pt_id'].unique())
	print(df.SV_type.unique())

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
	df = load_df('Z:/Members/clun/CLDB/meta/meta_LR.tsv')
	print(df)
	df_cluster = cluster(df)
	print('Writing Results to TSV')
	df_cluster.to_csv('CLDB_snf2_hg38.tsv', sep='\t', index=False)
	print('Done')
	
	### sqlite
	# fns = glob.glob('../data/CLDB_snf2/*hg19*anno*')
	# delete()
	# for fn in fns:
	# 	write_to_DB(fn)
	# create_index()
	# df = query()
	# print(df.columns)
	# print(df)
	# query()
if __name__ == '__main__':
	main()