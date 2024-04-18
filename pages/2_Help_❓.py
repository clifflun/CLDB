import streamlit as st

st.set_page_config(
    page_title="Help",
    page_icon="❓",
)

st.title("Q&A")

with st.expander("How many datasets are there?"):
    st.markdown("""
        This Variant Database operates under two separate datasets:
        [Carvalho Lab](https://pnri.org/carvalho-lab/) and [GREGoR Consortium](https://gregorconsortium.org/). 
        Each dataset curates and manages distinct cohort, catering to different research interests 
        and clinical applications. 
        """)
with st.expander("Which reference genome are the data aligned to?"):
    st.markdown("""
        Note that all data from Carvalho Lab are mapped to hg19, where data from GREGoR are mapped to hg38.
        """)

 