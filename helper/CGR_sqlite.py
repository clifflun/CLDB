import pandas as pd
import numpy as np
import sqlite3 as sq
from sqlalchemy import create_engine
import glob
import os
import re
from datetime import date



pd.set_option('display.expand_frame_repr', False)


db_name = 'GREGoR_052224'

def write_to_DB(fn):
	print('inserting')
	engine = create_engine('sqlite:///GREGoR_052224.sqlite')
	with pd.read_csv(fn, sep='\t', chunksize=50000, index_col=False) as reader:
		for chunk in reader:
			chunk.to_sql('CGR_hg38', con=engine, if_exists='append', index=False)
		

def create_index():
	conn = sq.connect('GREGoR_052224.sqlite')
	cursor = conn.cursor()
	_list = ['chr', 'type', 'cnv_lvl', 'pt_id', 'family', 'project', 'is_proband', 'cluster', 'SD_Overlap', 'OMIM_Count', 'match_CNV_SV', 'match_SV']
	for i in _list:
		print(i)
		cursor.execute(f'CREATE INDEX if NOT EXISTS {i}_index ON CGR_hg38 ({i})')
	conn.commit()
	conn.close()

def query():
	conn = sq.connect(f'GREGoR_052224.sqlite')
	df = pd.read_sql('select * from CGR_hg38 where 1=1 AND is_proband = 1 LIMIT 100 OFFSET 0', conn)
	conn.close()
	return df

def delete():
	conn = sq.connect('GREGoR_052224.sqlite')
	cursor = conn.cursor()	
	cursor.execute('DROP TABLE IF EXISTS CGR_hg38')
	conn.commit()
	conn.close()
	print('deleted table')


def update():
	conn = sq.connect('GREGoR_052224.sqlite')
	cursor = conn.cursor()
	
	cursor.execute('UPDATE CGR_hg38 SET is_proband = (CASE WHEN pt_id IN ("BH15692-2", "BH15692-3", "BH15693-2", "BH15694-2", "BH15694-3", "BH15695-2", "BH15695-3", "BH15696-2", "BH15696-3") THEN 0 ELSE is_proband END) ')
	conn.commit()
	conn.close()


def main():

	fns = glob.glob('./CGR_gregor/*')
	delete()
	for fn in fns:
		write_to_DB(fn)
	create_index()
	df = query()
	print(df.columns)
	print(df)



if __name__ == '__main__':
	main()