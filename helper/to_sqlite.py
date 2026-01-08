import pandas as pd
import numpy as np
import sqlite3 as sq
from sqlalchemy import create_engine, String
import glob
import os
import re
from datetime import date

pd.set_option('display.expand_frame_repr', False)

def write_meta_to_DB(fn):
	df = pd.read_csv(fn, sep='\t')
	df = df[['pt_id', 'sex', 'phenotype']]
	print(df)
	engine = create_engine('sqlite:///prod_CLDB_SR.sqlite')
	df.to_sql('meta_hg38', con=engine, if_exists='replace', index=False)
		
def write_to_DB():

	engine = create_engine('sqlite:///prod_CLDB_SR.sqlite')
	fnames = glob.glob('Z:/Members/clun/CLDB/CLDB_CNV_hg38_updated_annotated.tsv')
	print(fnames)
	for fn in fnames: 
		print(fn)
		with pd.read_csv(fn, sep='\t', chunksize=100000) as reader:
			for idx, chunk in enumerate(reader):
				print(f'processing chunk {idx}')
				chunk.to_sql('CNV_hg38', con=engine, if_exists='append', index=False)

def create_index():
	conn = sq.connect('prod_CLDB_SR.sqlite')
	cursor = conn.cursor()
	_list = ['chrom1', 'chrom2', 'SV_ID', 'SV_TYPE', 'PT_ID', 'FAM_ID', 'PROJECT', 'IS_PROBAND', \
	'CLUSTER_ID', 'COUNT', 'SD_overlap', 'OMIM_count', 'L_RefSeq', 'L_repeatmask', 'R_RefSeq', \
	'R_repeatmask', 'UNIQUE_PT_COUNT', 'imprinting_genes']
	for i in _list:
		print(i)
		cursor.execute(f'CREATE INDEX if NOT EXISTS {i}_index ON CNV_hg38 ({i})')
	conn.commit()
	conn.close()

def delete():
	print(f'deleting table')
	conn = sq.connect('prod_CLDB_SR.sqlite')
	cursor = conn.cursor()	
	cursor.execute('DROP TABLE IF EXISTS CNV_hg38')
	conn.commit()
	conn.close()
	print('deleted table')

def query():
	conn = sq.connect('prod_CLDB_SR.sqlite')
	print('running query')
	# df = pd.read_sql(f"""SELECT COUNT(*) from CNV_hg38 where 1=1 AND CLUSTER_ID = '63047' ORDER BY 1 LIMIT 100""", conn)
	df = pd.read_sql(f"""SELECT COUNT(*) from CNV_hg38 WHERE 1=1""", conn)
	df.to_csv('query_out3.csv', index=None)
	conn.close()
	return df

def update_column():
	df = pd.read_csv('Z:/Members/clun/CLDB/CLDB_CNV_hg38_updated_annotated.tsv', sep='\t')
	print(df.shape)
	conn = sq.connect('prod_CLDB_SR.sqlite')
	df.to_sql('gap_update', conn, if_exists='replace', index=False)
	print('updating')
	cursor = conn.cursor()
	update_query = """
	UPDATE CNV_hg38
	SET SD_overlap = (
	    SELECT SD_overlap
	    FROM gap_update
	    WHERE gap_update.UUID = CNV_hg38.UUID
	)
	WHERE UUID IN (SELECT uuid FROM gap_update)
	;
	"""

	cursor.execute(update_query)
	conn.commit()
	conn.close()



def main():
	delete()
	write_to_DB()
	create_index()
	write_meta_to_DB('./meta/meta_SR_hg38.tsv')
	update_column()
	df = query()
	print(df)


if __name__ == '__main__':
	main()
