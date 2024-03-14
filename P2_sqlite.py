import pandas as pd
import numpy as np
import sqlite3 as sq
from sqlalchemy import create_engine
import glob
import os
import re
from datetime import date



pd.set_option('display.expand_frame_repr', False)


db_name = 'P2_DB_070224'

def write_to_DB():

	engine = create_engine('sqlite:///test.sqlite')
	fnames = glob.glob('./070224/*annotated*')
	for fn in fnames: 
		print(fn)
		with pd.read_csv(fn, sep='\t', chunksize=50000) as reader:
			for chunk in reader:
				chunk.to_sql('P2_DB', con=engine, if_exists='append', index=False)
		

def create_index():
	conn = sq.connect('P2_DB_070224.sqlite')
	cursor = conn.cursor()
	_list = ['chrom1', 'chrom2', 'SV_id', 'SV_type', 'pt_id', 'family', 'project', 'is_proband', 'cluster', 'count', 'SD_Overlap', 'OMIM_Count', 'RefSeq_Disrupt_left', 'RefSeq_Disrupt_right', 'RepeatMask_Disrupt_left', 'RepeatMask_Disrupt_right', 'unique_pt_id_count']
	for i in _list:
		print(i)
		cursor.execute(f'CREATE INDEX if NOT EXISTS {i}_index ON P2_DB ({i})')
	conn.commit()
	conn.close()

def query():
	conn = sq.connect(f'P2_DB_070224.sqlite')
	df = pd.read_sql('select * from P2_DB where 1=1 AND is_proband = 1 LIMIT 100 OFFSET 0', conn)
	conn.close()
	return df


def update():
	conn = sq.connect('P2_DB_070224.sqlite')
	cursor = conn.cursor()
	
	cursor.execute('UPDATE P2_DB SET is_proband = (CASE WHEN pt_id IN ("BH15692-2", "BH15692-3", "BH15693-2", "BH15694-2", "BH15694-3", "BH15695-2", "BH15695-3", "BH15696-2", "BH15696-3") THEN 0 ELSE is_proband END) ')
	conn.commit()
	conn.close()


def query_multigene(gene_list):
	conn = sq.connect(f'{db_name}.sqlite')
	grouped = gene_list.groupby('inv_chr')
	for name, group in grouped: 
		if name != "X":
			continue
		dfs = []	
		genes = group.gene_symbol.unique()
		qry = 'SELECT * FROM P2_DB WHERE 1=1'
		qry += f' AND chrom1 = "{name}"'
		qry+= f' AND (RefSeq_Disrupt_left LIKE "%{genes[0]}%"'
		for gene in genes[1:]:
			qry += f' OR RefSeq_Disrupt_left LIKE "%{gene}%"'
		qry+= ')'
		print(qry)
		df = pd.read_sql(qry, conn)
		print(df.shape)
		
		qry = 'SELECT * FROM P2_DB WHERE 1=1'
		qry += f' AND chrom1 = "{name}"'
		qry+= f' AND (RefSeq_Disrupt_right LIKE "%{genes[0]}%"'
		for gene in genes[1:]:
			qry += f' OR RefSeq_Disrupt_right LIKE "%{gene}%"'
		qry+= ')'
		print(qry)
		df = pd.read_sql(qry, conn)
		print(df.shape)
		
		dfs.append(df)
		combined_df = pd.concat(dfs, ignore_index=True)
		combined_df_no_dup = combined_df.drop_duplicates(keep='first')

		combined_df_no_dup.to_csv(f'tugce_chr{name}.tsv', sep='\t', index=False)
	conn.close()


def main():
	# write_to_DB()
	create_index()
	# update()
	df = query()
	print(df.columns)
	print(df)

	# gene_list = pd.read_csv('Z:/Members/clun/analysis_PIDD/INV/DB_query/gene_list.txt', sep=' ') 
	# query_multigene(gene_list)
	


if __name__ == '__main__':
	main()