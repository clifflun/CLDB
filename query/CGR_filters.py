import streamlit as st
import pandas as pd
import numpy as np

def pheno2gene(hpo_id):
	hpo = pd.read_csv('./reference/phenotype_to_genes.txt', sep='\t')
	genes=[]
	out=hpo[hpo['hpo_id']==hpo_id]['gene_symbol'].unique().tolist()
	return out
	
chrom = [i for i in range(1,23)]
chrom.append('X')
chrom.append('Y')

OMIM = pd.read_csv('./reference/hg19/OMIM_gene2_hg19_UCSC_all.bed', sep='\t', index_col=False)
OMIM = OMIM[OMIM.pheno_key.notnull()]
OMIM_syms = OMIM.gene_symbol.unique()

rs = pd.read_csv(f'./reference/hg19/RefSeq_hg19.tsv', sep='\t', index_col=False)
rs_syms = rs.gene_id.unique()


def CGR_filters(dataset, key):
	if dataset == 'cldb':
		metadata = pd.read_csv('./meta/meta.tsv', sep='\t', index_col=False)
	elif dataset == 'gregor':
		metadata = pd.read_csv('./meta/gregor_meta.tsv', sep='\t', index_col=False)
	families = metadata.family.unique()
	pt_ids = metadata.pt_id.unique()
	projects = metadata.project.unique()
	n_sub = len(pt_ids)

	st.header('Currently querying:')
	st.write('Complex Genomic Rearrangement Calls')

	col1, col2, col3, col4, col5= st.columns(5)
	with col1:
		chrom1 = st.selectbox(
			'Chromosome',
			chrom,
			index=None,
			placeholder='Select Chromosome',
			key=f'chrom1-{key}'
			)

		project = st.multiselect(
		'Project:', 
		projects,
		key=f'project-{key}'
		)

		SD_ol = st.selectbox(
			'SD Overlap',
			['True', 'False'],
			index=None,
			help = 'If CGR call is overlapping >98% with Segmental Duplications (UCSC)',
			key=f'SD_ol-{key}'
			)

	with col2:
		start = st.number_input(
			'Left Breakpoint:', 
			value=None, 
			min_value=0,
			key=f'start-{key}')

		is_proband = st.selectbox(
		'Is call from a proband:',
		['True', 'False'],
		index=None,
		key=f'is_proband-{key}')

		cluster_id = st.number_input(
			'Cluster ID:', 
			value=None, 
			min_value=-1,
			key = f'cluster_id-{key}',
			help='CGR calls are clustered using DBSCAN (eps:500, min_sample:2)')

	with col3:
		end = st.number_input(
			'Right breakpoint:', 
			value=None, 
			min_value=0,
			key=f'end-{key}')

		pt_id = st.multiselect(
		'Individual ID:',
		pt_ids,
		key=f'pt_id-{key}')
	
		if st.checkbox('Filter CGR Length', key=f'svlen-{key}'):
			col_len1, col_len2 = st.columns(2)
			with col_len1:
				len_min = st.number_input(
					'min',
					min_value=0,
					key=f'svlen-min-{key}'
					)
			with col_len2:
				len_max = st.number_input(
					'max',
					min_value=len_min+1,
					key=f'svlen-max-{key}'
					)
		else: 
			len_min = None
			len_max = None
	with col4:
		svtype = st.selectbox(
			'CGR Type:',
			['gain', 'loss'],
			index=None,
			key=f'svtype-{key}'
			)

		family = st.multiselect(
		'Family ID:',
		families,
		help='Results will include all family members',
		key=f'family-{key}'
		)

	with col5:
		level = st.multiselect(
			'CGR level    ',
			['DUP', 'TRP', 'MUL_GAIN', 'HET_DEL', 'HOM_DEL', 'UND'],
			help='Characterized level of CGR call',
			key=f'level-{key}'
			)

		

	num_rows = st.selectbox(
	'Show Number of Rows:',
	[100,500,1000,5000,'All'],
	help='Use lower number of rows to test your query for faster query time',
	key=f'nrows-{key}')

	qry_dict={
	'chrom1': chrom1,
	'start': start,
	'end': end,
	'svtype': svtype,
	'level': level,
	'len_min': len_min, 
	'len_max': len_max,
	'SD_ol': SD_ol, 
	'is_proband': is_proband,
	'project': project, 
	'family': family,
	'pt_id': pt_id,
	'num_rows': num_rows,
	'cluster_id': cluster_id
	}
	return qry_dict


def CGR_build(qry_dict, table):
	
	qry = f'SELECT c.*, m.sex, m.phenotype FROM {table} c JOIN meta m ON c.pt_id = m.pt_id WHERE 1=1'

	if qry_dict['chrom1'] is not None:
				chrom1 = qry_dict['chrom1']
				qry += f' AND chr = "chr{chrom1}"'

	if qry_dict['start'] is not None:
		start = qry_dict['start']
		qry += f' AND start >= {start}'

	if qry_dict['end'] is not None:
		end = qry_dict['end']
		qry += f' AND end <= {end}'

	if len(qry_dict['level']) != 0:
		level = '", "'.join(qry_dict['level'])
		qry += f' AND cgr_lvl IN ("{level}")'

	if qry_dict['svtype'] is not None:
		svtype = qry_dict['svtype']
		qry += f' AND type = "{svtype}"'

	if qry_dict['len_min'] is not None:
		len_min = qry_dict['len_min']
		qry += f' AND len >= {len_min}'

	if qry_dict['len_max'] is not None:
		len_max = qry_dict['len_max']
		qry += f' AND len <= {len_max}'
	
	if qry_dict['SD_ol'] == 'True':
		SD_ol = 1
		qry += f' AND SD_Overlap = {SD_ol}'
	elif qry_dict['SD_ol'] == 'False':
		SD_ol = 0
		qry += f' AND SD_Overlap = {SD_ol}'

	if qry_dict['is_proband'] == 'True':
		is_proband = 1
		qry += f' AND is_proband = {is_proband}'
	elif qry_dict['is_proband'] == 'False':
		is_proband = 0
		qry += f' AND is_proband = {is_proband}'

	if len(qry_dict['project']) != 0:
		project = '", "'.join(qry_dict['project'])
		qry += f' AND project IN ("{project}")'

	if len(qry_dict['family']) != 0:
		family = '", "'.join(qry_dict['family'])
		qry += f' AND family IN ("{family}")'

	if len(qry_dict['pt_id']) != 0:
		pt_id = '", "'.join(qry_dict['pt_id'])
		qry += f' AND m.pt_id IN ("{pt_id}")'

	if qry_dict['cluster_id'] is not None:
		cluster_id = qry_dict['cluster_id']
		qry += f' AND cluster = {cluster_id}'

	qry += ' AND (left_dist <= 1500 OR right_dist <=1500) AND match_CNV_SV != 1'

	if qry_dict['num_rows'] == 'All':
		num_rows = 184467440737095516
	else:
		num_rows = qry_dict['num_rows']

	qry += f' LIMIT {num_rows}'
	return qry