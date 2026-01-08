import streamlit as st
from query.query import igv_query


st.set_page_config(
    page_title="IGV",
    page_icon="🧬",
    layout="wide"
)


sr_tab, lr_tab = st.tabs(['Short Read', 'Long Read'])
with sr_tab:
	sr_CC_hg19_tab, sr_CC_hg38_tab, sr_gg_hg38_tab = st.tabs(['Carvalho Lab(hg19)', 'Carvalho Lab(hg38)', 'GREGoR(hg38)'])
	with sr_CC_hg19_tab:
		igv_query('cldb', 'CNV', 'hg19', 'CLDB_SR')
	with sr_CC_hg38_tab:
		igv_query('cldb', 'CNV', 'hg38', 'CLDB_SR')
	with sr_gg_hg38_tab:
		igv_query('gregor', 'CNV', 'hg38', 'GREGoR_SR')

with lr_tab:
	lr_CC_hg19_tab, lr_CC_hg38_tab = st.tabs(['Carvalho Lab(hg19)', 'Carvalho Lab(hg38)'])
	with lr_CC_hg19_tab:
		igv_query('cldb', 'snf2', 'hg19', 'CLDB_LR')
	with lr_CC_hg38_tab:
		igv_query('cldb', 'snf2', 'hg38', 'CLDB_LR')




