import streamlit as st
from datetime import datetime

st.set_page_config(
    page_title="Hello",
    page_icon="👋",
)

print(f'Started CLDB at {datetime.now()}')

st.write("# Welcome to our Carvalho Lab Database (3.0.0)! 👋")

st.markdown(
"""
Welcome to our Database, a comprehensive repository of genetic variants gathered 
from various cohorts. Our database serves as a valuable resource for researchers, clinicians, 
and genetic counselors to explore and analyze genetic variations associated with diverse 
phenotypes and diseases.

**Getting Started:**

To begin exploring our Variant Database, simply set your filters to discover variants of interest. Whether you're investigating rare diseases, studying population genetics, or conducting pharmacogenomic research, our database provides valuable insights into the genetic basis of human traits and diseases.

We are committed to maintaining the accuracy, integrity, and usability of our database to support your research and clinical endeavors. If you have any questions, feedback, or suggestions, please don't hesitate to contact our team. Thank you for choosing our Variant Database as your trusted resource for genomic variant information.

Happy exploring!


**Having Issues?** 
 [Contact Us](https://github.com/clifflun/CLDB/issues)


"""
)
