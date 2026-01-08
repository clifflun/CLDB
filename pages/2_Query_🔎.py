import pandas as pd
import numpy as np
import streamlit as st
import plotly.express as px
from query.query import final_query

st.set_page_config(
    page_title="Query",
    page_icon="🔎",
    layout="wide"
)


sr_tab, lr_tab = st.tabs(['Short Read', 'Long Read'])
with sr_tab:
	sr_CC_hg19_tab, sr_CC_hg38_tab, sr_gg_hg38_tab = st.tabs(['Carvalho Lab(hg19)', 'Carvalho Lab(hg38)', 'GREGoR(hg38)'])
	with sr_CC_hg19_tab:
		# st.warning('Please be aware that sequencing quality is poor for BH16606-1, BH12526-1, BH12519-2. Calls might not be accurate.')
		sr_CC_hg19_tab_1, sr_CC_hg19_tab_2, sr_CC_hg19_tab_3 = st.tabs(['Parliament2', 'CNV', 'CGR'])
		with sr_CC_hg19_tab_1:
			final_query('cldb', 'P2', 'hg19', 'CLDB_SR')
		with sr_CC_hg19_tab_2:
			final_query('cldb', 'CNV', 'hg19', 'CLDB_SR')
		with sr_CC_hg19_tab_3:
			final_query('cldb', 'CGR', 'hg19', 'CLDB_SR')

	with sr_CC_hg38_tab:
		# st.warning('Please be aware that sequencing quality is poor for BH16606-1, BH12526-1, BH12519-2. Calls might not be accurate.')
		sr_CC_hg38_tab_1, sr_CC_hg38_tab_2, sr_CC_hg38_tab_3 = st.tabs(['Parliament2', 'CNV', 'CGR'])
		with sr_CC_hg38_tab_1:
			final_query('cldb', 'P2', 'hg38', 'CLDB_SR')
		with sr_CC_hg38_tab_2:
			final_query('cldb', 'CNV', 'hg38', 'CLDB_SR')
		with sr_CC_hg38_tab_3:
			final_query('cldb', 'CGR', 'hg38', 'CLDB_SR')

	with sr_gg_hg38_tab:
		sr_gg_hg38_tab_1, sr_gg_hg38_tab_2, sr_gg_hg38_tab_3 = st.tabs(['Parliament2', 'CNV', 'CGR'])
		with sr_gg_hg38_tab_1:			
			final_query('gregor', 'P2', 'hg38', 'GREGoR_SR')
		with sr_gg_hg38_tab_2:	
			final_query('gregor', 'CNV', 'hg38', 'GREGoR_SR')
		with sr_gg_hg38_tab_3:
			final_query('gregor', 'CGR', 'hg38', 'GREGoR_SR')

with lr_tab:
	lr_CC_hg19_tab, lr_CC_hg38_tab = st.tabs(['Carvalho Lab(hg19)', 'Carvalho Lab(hg38)'])
	with lr_CC_hg19_tab:
		final_query('cldb', 'snf2', 'hg19', 'CLDB_LR')
	with lr_CC_hg38_tab:
		final_query('cldb', 'snf2', 'hg38', 'CLDB_LR')

