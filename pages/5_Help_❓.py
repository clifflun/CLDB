import streamlit as st
import pandas as pd
import os

st.set_page_config(
	page_title="Help",
	page_icon="",
	)


st.title("Q&A")

df = pd.read_csv('prod_qa.txt', sep='\t', header=None)

def qa_generator(i,q,a):
	with st.expander(q):
		st.markdown(a)
		match i:
		    case 2:
		        st.image('./misc/P2_DB.png', caption="P2_SV pipeline")
		    case 4:
		        st.image('./misc/CGR.png', caption="A. Example of matching CNV and SV (filtered out). B. Example of non-matching CNV and SV (not filtered out)")


for row in df.itertuples():
	i=row.Index
	q=row[1]
	a=row[2]
	qa_generator(i, q, a)
