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
# st.dataframe(df)

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

cc, gregor = st.tabs(['Carvalho Lab', 'GREGoR'])
with cc:
	charts, tables = st.tabs(['Charts', 'Tables'])
	with charts:
		df = pd.read_csv('./meta/meta.tsv', sep='\t')
		demo_plots(df)
	with tables:
		st.dataframe(df.drop(['P2_path', 'MD_path'], axis=1))
with gregor:
	charts_gg, tables_gg = st.tabs(['Charts', 'Tables'])
	with charts_gg:
		df = pd.read_csv('./meta/gregor_meta.tsv', sep='\t')
		demo_plots(df)
	with tables_gg:
		st.dataframe(df.drop(['P2_path', 'MD_path'], axis=1))