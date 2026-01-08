from pyhpo import Ontology
import pandas as pd

# _ = Ontology()

# def get_term_name(term_id):
# 	term = Ontology.get_hpo_object(term_id)
# 	return term.name

# df = pd.read_csv('X:/Projects/GREGoR/meta/tmp.tsv', sep='\t')


# df['term_name']=df['term_id'].map(get_term_name)
# df.to_csv('X:/Projects/GREGoR/meta/phenotype_with_term_names.tsv', sep='\t', index=False)


df = pd.read_csv('W:/Projects/GREGoR/meta/phenotype_with_term_names.tsv', sep='\t')

df2=df.pivot_table(index="participant_id", values=["term_name", "term_id"], aggfunc=lambda x: ", ".join(x))

df2=df2.reset_index()

df2.to_csv('W:/Projects/GREGoR/meta/phenotype_with_term_names_agg.tsv', sep='\t', index=False)
