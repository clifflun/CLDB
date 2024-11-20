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


def CNV_filters(dataset, key):
	if dataset == 'cldb':
		metadata = pd.read_csv('./meta/meta.tsv', sep='\t', index_col=False)
	elif dataset == 'gregor':
		metadata = pd.read_csv('./meta/gregor_meta.tsv', sep='\t', index_col=False)
	families = metadata.family.unique()
	pt_ids = metadata.pt_id.unique()
	projects = metadata.project.unique()
	n_sub = len(pt_ids)

	st.header('Currently querying:')
	st.write('Read depth-based CNV Calls')
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
			help = 'If CNV call is overlapping >98% with Segmental Duplications (UCSC)',
			key=f'SD_ol-{key}'
			)
		if st.checkbox('Filter CNV Length', key=f'svlen-{key}'):
			col_sv_len1, col_sv_len2 = st.columns(2)
			with col_sv_len1:
				sv_len_min = st.number_input(
					'min',
					min_value=0,
					key=f'svlen-min-{key}'
					)
			with col_sv_len2:
				sv_len_max = st.number_input(
					'max',
					min_value=sv_len_min+1,
					key=f'svlen-max-{key}'
					)
		else: 
			sv_len_min = None
			sv_len_max = None

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
			help='CNV calls are clustered using DBSCAN (eps:500, min_sample:2)')


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


		hpo_list = st.text_input(
			'HPO terms Search', '',
			key=f'hpo_genes-{key}',
			help='comma-separated HPO terms list (e.g. "HP:0003300, HP:0000200")'
			)

		

	with col4:
		cnvtype = st.selectbox(
			'CNV Type:',
			['gain', 'loss'],
			index=None,
			key=f'cnvtype-{key}'
			)

		family = st.multiselect(
		'Family ID:',
		families,
		help='Results will include all family members',
		key=f'family-{key}'
		)

		gene_list = st.text_input(
			'RefSeq Gene List Search', '',
			key=f'gene_list-{key}',
			help='comma-separated gene list (e.g. "MECP2, DOCK8")'
			)

	with col5:
		level = st.multiselect(
			'CNV level    ',
			['DUP', 'TRP', 'MUL_GAIN', 'HET_DEL', 'HOM_DEL', 'UND'],
			help='Characterized level of CNV call',
			key=f'level-{key}'
			)

		OMIM_sym = st.selectbox(
		'OMIM Disease Genes:',
		OMIM_syms,
		index=None,
		help='OMIM disease causing genes from UCSC (pheno_key 3/4)',
		key=f'OMIMsyms-{key}'
		)

		rs_sym = st.selectbox(
		'RefSeq Genes:',
		rs_syms,
		index=None,
		help='RefSeq genes from NCBI curated list',
		key=f'rssym-{key}'
		)


	col6, col7 = st.columns(2)
	with col6: 
		with st.expander('Breakpoint Filters'):	


			if st.checkbox('Left Breakpoint Disrupting Gene',
				help='check if left breakpoint is distrupting NCBI curated gene list',
				key=f'left-bkpt-gene-{key}'):
				disrupt_gene_left = True
			else: 
				disrupt_gene_left = False

			if st.checkbox('Right Breakpoint Disrupting Gene',
				help='check if right breakpoint is distrupting NCBI curated gene list',
				key=f'right-bkpt-gene-{key}'):
				disrupt_gene_right = True
			else: 
				disrupt_gene_right = False


			if st.checkbox('Left Breakpoint In Repeats',
				help='check if left breakpoint is distrupting RepeatMask repeats',
				key=f'left-bkpt-repeat-{key}'):
				disrupt_repeat_left = True
			else: 
				disrupt_repeat_left = False

			if st.checkbox('Right Breakpoint In Repeats',
				help='check if right breakpoint is distrupting RepeatMask repeats',
				key=f'right-bkpt-repeat-{key}'):
				disrupt_repeat_right = True
			else: 
				disrupt_repeat_right = False

		with st.expander('Dosage Sensitivity Filters'):
			if st.checkbox('CNV encompasses haploinsufficient genes (Collins)  ',
				help='HI genes from Collins et al 2022',
				key=f'hi-collins-{key}'):
				pHaplo_collins = True
			else: 
				pHaplo_collins = False

			if st.checkbox('CNV encompasses triplosensitive genes (Collins) ',
				help='TS genes from Collins et al 2022',
				key=f'ts-collins-{key}'):
				pTriplo_collins = True
			else: 
				pTriplo_collins = False

			if st.checkbox('CNV encompasses haploinsufficient genes (ClinGen) ',
				help='HI genes from ClinGen(2024)',
				key=f'hi-clingen-{key}'):
				pHaplo_clingen = True
			else: 
				pHaplo_clingen = False

			if st.checkbox('CNV encompasses triplosensitive genes (ClinGen) ',
				help='TS genes from ClinGen(2024)',
				key=f'ts-clingen-{key}'):
				pTriplo_clingen = True
			else: 
				pTriplo_clingen = False

			if st.checkbox('CNV encompasses DECIPHER CNV Syndromes regions ',
				help='Syndromes list from DECIPHER',
				key=f'decipher-{key}'):
				decipher_cnv = True
			else: 
				decipher_cnv = False

			if st.checkbox('CNV overlaps ISCA regions    ',
				help='ISCA regions',
				key=f'isca-{key}'):
				ISCA_cnv = True
			else: 
				ISCA_cnv = False
	with col7: 
		with st.expander('Count/Frequency Filters'):	

			if st.checkbox('Filter pseudo-Database Frequency',
				help='Frequency of number of unique individuals in a cluster in the whole DB',
				key=f'db-freq-{key}'):
										
				db_freq = st.slider(
					' ',
					value=(0.0, 1.0),
					max_value=1.0,
					step=0.0001,
					format='%f',
					key=f'df-freq-slider-{key}')
			else: 
				db_freq = None

			if st.checkbox('Filter Unique Individual Count', 
				help='Number of unique individuals in a cluster in the whole DB',
				key=f'df-count-{key}'):
				count = st.slider(
					'        ',
					value=(1,n_sub),
					max_value=n_sub,
					key=f'df-count-slider-{key}'
					)
			else: 
				count = None

			if st.checkbox('Filter gnomAD Frequency',
				help='gnomAD v2 population frequency (if matched)',
				key=f'gnomad-freq-{key}'):
				gnomad_freq_value=(0.0, 1.0)
				gnomAD_freq = st.slider(
					'',
					value=gnomad_freq_value,
					max_value=1.0,
					step=0.0001,
					format='%f',
					key=f'gnomad-freq-slider-{key}')
			else: 
				gnomAD_freq = None	
			
			if st.checkbox('Filter RefSeq Count',
				help='Number of NCBI curated RefSeq genes overlapping a call',
				key=f'refseq-count-{key}'):
				col_rs1, col_rs2 = st.columns(2)
				with col_rs1:
					RefSeq_min = st.number_input(
						'min',
						min_value=0,
						key=f'refseq-count-min-{key}'
						)
				with col_rs2:
					RefSeq_max = st.number_input(
						'max',
						min_value=RefSeq_min+1,
						key=f'refseq-count-max-{key}'
						)
			else: 
				RefSeq_min = None
				RefSeq_max = None

			if st.checkbox('Filter OMIM Count',
				help='Number of OMIM disease-causing genes overlapping a call',
				key=f'OMIM-count-{key}'):
				col_omim1, col_omim2 = st.columns(2)
				with col_omim1:
					OMIM_min = st.number_input(
						'min',
						min_value=0,
						key=f'OMIM-count-min-{key}'
						)
				with col_omim2:
					OMIM_max = st.number_input(
						'max',
						min_value=OMIM_min+1,
						key=f'OMIM-count-max-{key}'
						)
			else: 
				OMIM_min = None
				OMIM_max = None


			

	num_rows = st.selectbox(
	'Show Number of Rows:',
	[100,500,1000,5000,'All'],
	help='Use lower number of rows to test your query for faster query time',
	key=f'nrows-{key}')

	qry_dict={
		'chrom1': chrom1,
		'start': start,
		'end': end,
		'cnvtype': cnvtype,
		'level': level,
		'sv_len_min': sv_len_min, 
		'sv_len_max': sv_len_max,
		'SD_ol': SD_ol, 
		'is_proband': is_proband,
		'count': count,
		'db_freq': db_freq,
		'gnomAD_freq': gnomAD_freq,
		'project': project, 
		'OMIM_min': OMIM_min,
		'OMIM_max': OMIM_max,
		'RefSeq_min': RefSeq_min,
		'RefSeq_max': RefSeq_max,
		'OMIM_sym': OMIM_sym,
		'rs_sym': rs_sym,
		'rs_list': gene_list,
		'hpo_list': hpo_list,
		'family': family,
		'pt_id': pt_id,
		'num_rows': num_rows,
		'cluster_id': cluster_id,
		'disrupt_gene_left': disrupt_gene_left,
		'disrupt_gene_right': disrupt_gene_right,
		'disrupt_repeat_left': disrupt_repeat_left,
		'disrupt_repeat_right': disrupt_repeat_right,
		'pHaplo_collins': pHaplo_collins,
		'pTriplo_collins': pTriplo_collins,
		'pHaplo_clingen': pHaplo_clingen,
		'pTriplo_clingen': pTriplo_clingen,
		'decipher_cnv': decipher_cnv,
		'ISCA_cnv': ISCA_cnv}
	return qry_dict


def CNV_build(qry_dict, table):
	if len(qry_dict['hpo_list']) != 0:
		tmp = []
		for hpo in [x.strip() for x in qry_dict['hpo_list'].split(',')]:
			tmp.append(pheno2gene(hpo))
		out=list(set([x for xs in tmp for x in xs]))
		qry_dict['rs_list']=','.join(out)
	
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
		qry += f' AND cnv_lvl IN ("{level}")'

	if qry_dict['cnvtype'] is not None:
		cnvtype = qry_dict['cnvtype']
		qry += f' AND type = "{cnvtype}"'

	if qry_dict['sv_len_min'] is not None:
		sv_len_min = qry_dict['sv_len_min']
		qry += f' AND len >= {sv_len_min}'

	if qry_dict['sv_len_max'] is not None:
		sv_len_max = qry_dict['sv_len_max']
		qry += f' AND len <= {sv_len_max}'
	
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

	if qry_dict['count'] is not None:
		count1 = qry_dict['count'][0]
		count2 = qry_dict['count'][1]	
		qry += f' AND unique_pt_id_count BETWEEN {count1} and {count2}'

	if qry_dict['db_freq'] is not None:
		pf1 = qry_dict['db_freq'][0]
		pf2 = qry_dict['db_freq'][1]	
		qry += f' AND psuedo_df_freq BETWEEN {pf1} and {pf2}'

	if qry_dict['gnomAD_freq'] is not None:
		gf1 = qry_dict['gnomAD_freq'][0]
		gf2 = qry_dict['gnomAD_freq'][1]	
		qry += f' AND gnomAD_AF BETWEEN {gf1} and {gf2}'

	if len(qry_dict['project']) != 0:
		project = '", "'.join(qry_dict['project'])
		qry += f' AND project IN ("{project}")'

	if qry_dict['OMIM_min'] is not None:
		OMIM_min = qry_dict['OMIM_min']
		qry += f' AND OMIM_Count >= {OMIM_min}'

	if qry_dict['OMIM_max'] is not None:
		OMIM_max = qry_dict['OMIM_max']
		qry += f' AND OMIM_Count <= {OMIM_max}'


	if qry_dict['RefSeq_min'] is not None:
		RefSeq_min = qry_dict['RefSeq_min']
		qry += f' AND RefSeq_Count >= {RefSeq_min}'

	if qry_dict['RefSeq_max'] is not None:
		RefSeq_max = qry_dict['RefSeq_max']
		qry += f' AND RefSeq_Count <= {RefSeq_max}'
	
	if qry_dict['OMIM_sym'] is not None:
		OMIM_sym = qry_dict['OMIM_sym']
		qry += f' AND OMIM_Symbol LIKE "%{OMIM_sym}%"'

	if qry_dict['rs_sym'] is not None:
		rs_sym = qry_dict['rs_sym']
		qry += f' AND RefSeq_Symbol LIKE "%{rs_sym}%"'


	if len(qry_dict['rs_list']) != 0:
		genes=[x.strip().upper() for x in qry_dict['rs_list'].split(',')]
		genes_2 = [genes[i:i+900] for i in range(0, len(genes), 900)]
		# st.write(genes_2)

		qry+=' AND ('
		for gene in genes_2:
			qry+= f'(RefSeq_Symbol LIKE "%{gene[0]}%"'
			for g in gene[1:]:
				qry += f' OR RefSeq_Symbol LIKE "%{g}%"'
			qry+= ') OR '
		qry+='RefSeq_Symbol LIKE "")'


	if len(qry_dict['family']) != 0:
		family = '", "'.join(qry_dict['family'])
		qry += f' AND family IN ("{family}")'

	if len(qry_dict['pt_id']) != 0:
		pt_id = '", "'.join(qry_dict['pt_id'])
		qry += f' AND m.pt_id IN ("{pt_id}")'


	if qry_dict['cluster_id'] is not None:
		cluster_id = qry_dict['cluster_id']
		qry += f' AND cluster = {cluster_id}'

	bkpt_filters = [qry_dict['disrupt_gene_left'], qry_dict['disrupt_gene_right'], qry_dict['disrupt_repeat_left'], qry_dict['disrupt_repeat_right']]
	if True in bkpt_filters:
		tmp_qry=[]
		if qry_dict['disrupt_gene_left']: 
			tmp_qry.append('RefSeq_Disrupt_left IS NOT NULL')

		if qry_dict['disrupt_gene_right']: 
			tmp_qry.append('RefSeq_Disrupt_right IS NOT NULL')

		if qry_dict['disrupt_repeat_left']: 
			tmp_qry.append('RepeatMask_Disrupt_left IS NOT NULL')

		if qry_dict['disrupt_repeat_right']: 
			tmp_qry.append('RepeatMask_Disrupt_right IS NOT NULL')

		qry += f' AND ({" OR ".join(tmp_qry)})'

	dosage_filters = [qry_dict['pHaplo_collins'], qry_dict['pTriplo_collins'], qry_dict['pHaplo_clingen'], qry_dict['pTriplo_clingen'], qry_dict['decipher_cnv'], qry_dict['ISCA_cnv']]
	if True in dosage_filters:
		tmp_qry = []

		if qry_dict['pHaplo_collins']: 
			tmp_qry.append('pHaplo_Collins IS NOT NULL')

		if qry_dict['pTriplo_collins']: 
			tmp_qry.append('pTriplo_Collins IS NOT NULL')

		if qry_dict['pHaplo_clingen']: 
			tmp_qry.append('pHaplo_clinGen IS NOT NULL')

		if qry_dict['pTriplo_clingen']: 
			tmp_qry.append('pTriplo_clinGen IS NOT NULL')

		if qry_dict['decipher_cnv']: 
			tmp_qry.append('DECIPHER_CNV_Syndromes IS NOT NULL')
		if qry_dict['ISCA_cnv']: 	
			tmp_qry.append('ISCA IS NOT NULL')

		qry += f' AND ({" OR ".join(tmp_qry)})'


	if qry_dict['num_rows'] == 'All':
		num_rows = 184467440737095516
	else:
		num_rows = qry_dict['num_rows']
	qry += f' LIMIT {num_rows}'
	return qry


