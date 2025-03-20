import pandas as pd
import numpy as np
import sqlite3 as sq
from sqlalchemy import create_engine
import glob
import os
import re
from datetime import date

pd.set_option('display.expand_frame_repr', False)

def write_meta_to_DB(fn):
	df = pd.read_csv(fn, sep='\t')
	df = df[['pt_id', 'sex', 'phenotype']]
	print(df)
	engine = create_engine('sqlite:///prod_CLDB_LR.sqlite')
	df.to_sql('meta', con=engine, if_exists='append', index=False)
		
def write_to_DB():

	engine = create_engine('sqlite:///prod_CLDB_LR.sqlite')
	fnames = glob.glob('Z:/Members/clun/CLDB/data/3.0.0/CLDB_snf2_hg19_annotated.tsv')
	# print(fnames)
	for fn in fnames: 
		print(fn)
		with pd.read_csv(fn, sep='\t', chunksize=100000) as reader:
			for chunk in reader:
				chunk.to_sql('snf2_hg19', con=engine, if_exists='append', index=False)

def create_index():
	conn = sq.connect('prod_CLDB_LR.sqlite')
	cursor = conn.cursor()
	_list = ['chrom1', 'chrom2', 'SV_ID', 'SV_TYPE', 'PT_ID', 'FAM_ID', 'PROJECT', 'IS_PROBAND', 'CLUSTER_ID', 'COUNT', 'SD_overlap', 'OMIM_count', 'L_RefSeq', 'L_repeatmask', 'R_RefSeq', 'R_repeatmask', 'UNIQUE_PT_COUNT', 'imprinting_genes']
	for i in _list:
		print(i)
		cursor.execute(f'CREATE INDEX if NOT EXISTS {i}_index ON snf2_hg19 ({i})')
	conn.commit()
	conn.close()

def delete():
	conn = sq.connect('prod_CLDB_LR.sqlite')
	cursor = conn.cursor()	
	cursor.execute('DROP TABLE IF EXISTS snf2_hg19')
	conn.commit()
	conn.close()
	print('deleted table')

def query():
	conn = sq.connect('prod_CLDB_LR.sqlite')
	print('running query')
	df = pd.read_sql('select * from snf2_hg19 where 1=1 ', conn)

	conn.close()
	return df

def update():
	print('updating')
	conn = sq.connect('prod_CLDB_LR.sqlite')
	cursor = conn.cursor()
	# cursor.execute('UPDATE snf2_hg19 SET project = (CASE WHEN pt_id IN ("BH13842-1", "BH13842-2", "BH13842-3") THEN "MECP2(-)Rett" ELSE project END) ')
	conn.commit()
	conn.close()

def main():
	delete()
	write_to_DB()
	create_index()
	write_meta_to_DB('./meta/meta_LR.tsv')
	# update()
	# df = query()
	# print(df)

if __name__ == '__main__':
	main()