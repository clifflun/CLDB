import pandas as pd
import numpy as np
import sqlite3 as sq
from sqlalchemy import create_engine
import glob
import os
import re
from datetime import date



pd.set_option('display.expand_frame_repr', False)


db_name = 'prod_GREGoR.sqlite'

def write_to_DB(fn):
	df = pd.read_csv(fn, sep='\t')
	df = df[['pt_id', 'sex', 'phenotype']]
	print(df)
	engine = create_engine('sqlite:///prod_GREGoR.sqlite')
	df.to_sql('GG_meta', con=engine, if_exists='append', index=False)
		

def create_index():
	conn = sq.connect('prod_GREGoR.sqlite')
	cursor = conn.cursor()
	_list = ['pt_id', 'sex', 'phenotype']
	for i in _list:
		print(i)
		cursor.execute(f'CREATE INDEX if NOT EXISTS {i}_index ON GG_meta ({i})')
	conn.commit()
	conn.close()

def query():
	conn = sq.connect('../prod_CLDB.sqlite')
	df = pd.read_sql('select * from meta where 1=1 LIMIT 100 OFFSET 0', conn)
	conn.close()
	return df

def delete():
	conn = sq.connect('prod_GREGoR.sqlite')
	cursor = conn.cursor()	
	cursor.execute('DROP TABLE IF EXISTS GG_meta')
	conn.commit()
	conn.close()
	print('deleted table')


def update():
	conn = sq.connect('../prod_CLDB.sqlite')
	cursor = conn.cursor()
	cursor.execute('ALTER TABLE "CC_meta" RENAME TO "meta"')
	conn.commit()
	conn.close()


def main():

	# delete()
	# write_to_DB('./meta/gregor_meta.tsv')
	# create_index()
	# df=query()
	# print(df)
	update()

if __name__ == '__main__':
	main()