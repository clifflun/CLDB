import pandas as pd
import numpy as np
import sqlite3 as sq
from sqlalchemy import create_engine
import glob
import os
import re
from datetime import date



pd.set_option('display.expand_frame_repr', False)


db_name = 'P2_DB_040324'

def write_to_DB(fn):

	engine = create_engine('sqlite:///P2_DB_040324.sqlite')
	with pd.read_csv(fn, sep='\t', chunksize=50000, index_col=False) as reader:
		for chunk in reader:
			chunk.to_sql('CNV_hg19', con=engine, if_exists='append', index=False)
		

def create_index():
	conn = sq.connect('P2_DB_040324.sqlite')
	cursor = conn.cursor()
	_list = ['chr', 'type', 'cnv_lvl', 'pt_id', 'family', 'project', 'is_proband', 'cluster', 'count', 'SD_Overlap', 'OMIM_Count', 'RefSeq_Disrupt_left', 'RefSeq_Disrupt_right', 'RepeatMask_Disrupt_left', 'RepeatMask_Disrupt_right', 'unique_pt_id_count']
	for i in _list:
		print(i)
		cursor.execute(f'CREATE INDEX if NOT EXISTS {i}_index ON CNV_hg19 ({i})')
	conn.commit()
	conn.close()

def query():
	conn = sq.connect(f'P2_DB_040324.sqlite')
	df = pd.read_sql('select * from CNV_hg19 where 1=1 AND is_proband = 1 LIMIT 100 OFFSET 0', conn)
	conn.close()
	return df

def delete():
	conn = sq.connect('P2_DB_040324.sqlite')
	cursor = conn.cursor()	
	cursor.execute('DROP TABLE IF EXISTS CNV_hg19')
	conn.commit()
	conn.close()
	print('deleted table')


def update():
	conn = sq.connect('P2_DB_040324.sqlite')
	cursor = conn.cursor()
	
	cursor.execute('UPDATE CNV_hg19 SET is_proband = (CASE WHEN pt_id IN ("BH15692-2", "BH15692-3", "BH15693-2", "BH15694-2", "BH15694-3", "BH15695-2", "BH15695-3", "BH15696-2", "BH15696-3") THEN 0 ELSE is_proband END) ')
	conn.commit()
	conn.close()


def main():
	delete()
	write_to_DB('cnv_DB_040424.tsv')
	create_index()
	df = query()
	print(df.columns)
	print(df)



if __name__ == '__main__':
	main()