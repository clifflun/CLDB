import pandas as pd

df = pd.read_csv('NCBI_RefSeq_hg19_clean.bed', sep='\t')
out = df.groupby(['gene_id', 'seqname']).agg({'start': 'min', 'end': 'max'}).reset_index()



out_reord = out[['seqname', 'start', 'end', 'gene_id']]

out_reord.rename(columns={out_reord.columns[0]: 'chrom'}, inplace=True)
out_reord['chrom'] = out_reord['chrom'].str.replace('chr','')


out_reord.to_csv('RefSeq_hg19.tsv', sep='\t', index=False)