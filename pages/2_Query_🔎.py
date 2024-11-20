import pandas as pd
import numpy as np
import streamlit as st
import plotly.express as px
from query.P2_filters import P2_filters, P2_build
from query.CNV_filters import CNV_filters, CNV_build
from query.CGR_filters import CGR_filters, CGR_build

st.set_page_config(
    page_title="Query",
    page_icon="🔎",
    layout="wide"
)

def query(qry, db):
	conn = st.connection(db, type='sql')
	df=conn.query(qry)
	return df   

def final_query(dataset, pipeline, ref, db):
	key=dataset+'_'+pipeline+'_'+ref
	table=pipeline+'_'+ref
	if pipeline == "P2":
		qry_dict=P2_filters(dataset, key)
		qry=P2_build(qry_dict, table)
	elif pipeline == "CNV":
		qry_dict=CNV_filters(dataset, key)
		qry=CNV_build(qry_dict, table)
	elif pipeline == "CGR":
		qry_dict=CGR_filters(dataset, key)
		qry=CGR_build(qry_dict, table)
	if st.button('Query', type='primary', key=f'qry-btn-{key}'):
		df=query(qry, db)
		st.markdown("""---""")
		st.write(f'{len(df.index)} results')
		st.dataframe(
			df,
			hide_index=False,
			height=600)	

CC_tab, gregor_tab = st.tabs(['Carvalho Lab(hg19)', 'GREGoR(hg38)'])

with CC_tab:
	st.warning('Please be aware that sequencing quality is poor for BH16606-1, BH12526-1, BH12519-2. Calls might not be accurate.')
	CC_tab1, CC_tab2, CC_tab3 = st.tabs(['Parliament2', 'CNV', 'CGR'])
	with CC_tab1:
		final_query('cldb', 'P2', 'hg19', 'CLDB')
	with CC_tab2:
		final_query('cldb', 'CNV', 'hg19', 'CLDB')
	with CC_tab3:
		final_query('cldb', 'CGR', 'hg19', 'CLDB')

with gregor_tab:
	gregor_tab1, gregor_tab2, gregor_tab3 = st.tabs(['Parliament2', 'CNV', 'CGR'])
	with gregor_tab1:			
		final_query('gregor', 'P2', 'hg38', 'GREGoR')
	with gregor_tab2:	
		final_query('gregor', 'CNV', 'hg38', 'GREGoR')
	with gregor_tab3:
		final_query('gregor', 'CGR', 'hg38', 'GREGoR')
		