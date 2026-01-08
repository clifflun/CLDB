import pandas as pd
import numpy as np
import streamlit as st
import plotly.express as px


st.set_page_config(
	page_title="Datasets",
	page_icon="🗃️",
	layout="wide"
)

st.title('Dataset Summary')



def demo_plots(df):
	col1, col2 = st.columns(2)
	with col1:
		projs = df.project.value_counts()
		fig = px.pie(names=projs.index, values=projs.values,title='Individual per project')
		fig.update_traces(textposition='inside')
		fig.update_layout(uniformtext_minsize=12, uniformtext_mode='hide')
		st.plotly_chart(fig, use_container_width=True)
	with col2:
		sex=df.groupby('project')['sex'].value_counts().to_frame().reset_index()
		fig = px.bar(sex, x='count', y='project', color='sex', title='Sex distribution')
		st.plotly_chart(fig, use_container_width=True)

	pheno = df.phenotype.value_counts()
	fig = px.pie(
		names=pheno.index, 
		values=pheno.values,
		color_discrete_sequence=px.colors.sequential.RdBu,
		title='Phenotypes')
	fig.update_traces(textposition='inside')
	fig.update_layout(uniformtext_minsize=12, uniformtext_mode='hide')
	st.plotly_chart(fig, use_container_width=True)

	is_proband=df.groupby('project')['is_proband'].value_counts().to_frame().reset_index()
	fig = px.bar(is_proband, x='count', y='project', color='is_proband', title='Proband')
	st.plotly_chart(fig, use_container_width=True)

def charts_tables(df, last_col):
	charts, tables = st.tabs(['Charts', 'Tables'])
	with charts:
		demo_plots(df)
	with tables:
		st.dataframe(df.loc[:,:f'{last_col}'])


cc_sr, cc_lr, gregor = st.tabs(['Carvalho Lab (SR)', 'Carvalho Lab (LR)', 'GREGoR (SR)'])
with cc_sr:
	cc_sr_hg19, cc_sr_hg38 = st.tabs(['hg19', 'hg38'])
	with cc_sr_hg19:
		df = pd.read_csv('/CLDB/meta/meta_SR_hg19.tsv', sep='\t')
		charts_tables(df,'Inheritance')	
	with cc_sr_hg38:
		st.write('Coming Soon')
with cc_lr:
	cc_lr_hg19, cc_lr_hg38 = st.tabs(['hg19', 'hg38'])
	with cc_lr_hg19:
		df = pd.read_csv('/CLDB/meta/meta_LR.tsv', sep='\t')
		df=df.loc[df['has_LR_hg19']==1]
		df.reset_index(inplace=True)
		charts_tables(df,'Inheritance')	
	with cc_lr_hg38:
		df = pd.read_csv('/CLDB/meta/meta_LR.tsv', sep='\t')
		df=df.loc[df['has_LR_hg38']==1]
		df.reset_index(inplace=True)
		charts_tables(df,'Inheritance')	
with gregor:
	gregor_hg19, gregor_hg38 = st.tabs(['hg19', 'hg38'])
	with gregor_hg19:
		st.write('Coming Soon')
	with gregor_hg38:
		df = pd.read_csv('/CLDB/meta/gregor_meta.tsv', sep='\t')
		charts_tables(df,'Pre-discovery OMIM disorders')	