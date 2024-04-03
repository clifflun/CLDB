import pandas as pd
import numpy as np
import streamlit as st
import plotly.express as px

st.set_page_config(layout="wide")
st.title('Carvalho Lab Variant Database (v1.2.0)')

chrom = [i for i in range(1,23)]
chrom.append('X')
chrom.append('Y')

metadata = pd.read_csv('./meta/P2_meta.tsv', sep='\t', index_col=False)
families = metadata.family.unique()
pt_ids = metadata.pt_id.unique()
projects = metadata.project.unique()
n_sub = len(pt_ids)

OMIM = pd.read_csv('./reference/OMIM_gene2_hg19_UCSC_all.bed', sep='\t', index_col=False)
OMIM = OMIM[OMIM.pheno_key.notnull()]
OMIM_syms = OMIM.gene_symbol.unique()

rs = pd.read_csv(f'./reference/RefSeq_hg19.tsv', sep='\t', index_col=False)
rs_syms = rs.gene_id.unique()

CC_tab, gregor_tab = st.tabs(['Cavalho Lab', 'GREGoR'])

with CC_tab:
	CC_tab1, CC_tab2, CC_tab3 = st.tabs(['P2_SV', 'CNV', 'CGR'])

	with CC_tab1:
		with st.container():
			st.header('Currently querying:')
			st.write('Parliament 2 SV Calls (hg19)')
			col1, col2, col3, col4, col5= st.columns(5)

			with col1:
				chrom1 = st.selectbox(
					'Chromosome',
					chrom,
					index=None,
					placeholder='Select Chromosome'
					)

				project = st.multiselect(
				'Project:', 
				projects,
				key = 'project'
				)

				SD_ol = st.selectbox(
					'SD Overlap',
					['True', 'False'],
					index=None,
					help = 'If SV call is overlapping >98% with Segmental Duplications (UCSC)')


			with col2:
				start = st.number_input(
					'Left Breakpoint:', 
					value=None, 
					min_value=0,
					key='start')

				is_proband = st.selectbox(
				'Is SV call from a proband:',
				['True', 'False'],
				index=None,
				key = 'is_proband')

				cluster_id = st.number_input(
					'Cluster ID:', 
					value=None, 
					min_value=-1,
					key = 'cluster_id',
					help='SV calls are clustered using DBSCAN (eps:500, min_sample:2)')


			with col3:
				end = st.number_input(
					'Right breakpoint:', 
					value=None, 
					min_value=0,
					key='end')

				pt_id = st.multiselect(
				'Individual ID:',
				pt_ids,
				key='pt_id')


				if st.checkbox('Filter SV Length'):
					col_sv_len1, col_sv_len2 = st.columns(2)
					with col_sv_len1:
						sv_len_min = st.number_input(
							'  min',
							min_value=0
							)
					with col_sv_len2:
						sv_len_max = st.number_input(
							'  max',
							min_value=sv_len_min+1
							)
				else: 
					sv_len_min = None
					sv_len_max = None
				

			with col4:
				svtype = st.multiselect(
					'SV Type:',
					['DUP', 'DEL', 'INV', 'INS', 'BND']
					)

				family = st.multiselect(
				'Family ID:',
				families,
				help='Results will include all family members')




			with col5:
				genotype = st.selectbox(
					'Genotype of SV call:', 
					['Homozygous', 'Heterozygous'],
					index=None)

				OMIM_sym = st.selectbox(
				'OMIM Disease Genes:',
				OMIM_syms,
				index=None,
				help='OMIM disease causing genes from UCSC (pheno_key 3/4)'
				)

				rs_sym = st.selectbox(
				'RefSeq Genes:',
				rs_syms,
				index=None,
				help='RefSeq genes from NCBI curated list'
				)


			col6, col7 = st.columns(2)
			with col6: 
				with st.expander('Breakpoint Filters'):	


					if st.checkbox('Left Breakpoint Disrupting Gene',
						help='check if left breakpoint is distrupting NCBI curated gene list'):
						disrupt_gene_left = True
					else: 
						disrupt_gene_left = False

					if st.checkbox('Right Breakpoint Disrupting Gene',
						help='check if right breakpoint is distrupting NCBI curated gene list'):
						disrupt_gene_right = True
					else: 
						disrupt_gene_right = False


					if st.checkbox('Left Breakpoint In Repeats',
						help='check if left breakpoint is distrupting RepeatMask repeats'):
						disrupt_repeat_left = True
					else: 
						disrupt_repeat_left = False

					if st.checkbox('Right Breakpoint In Repeats',
						help='check if right breakpoint is distrupting RepeatMask repeats'):
						disrupt_repeat_right = True
					else: 
						disrupt_repeat_right = False

			with col7: 
				with st.expander('Count/Frequency Filters'):	

					if st.checkbox('Filter pseudo-Database Frequency',
						help='Frequency of number of unique individuals in a cluster in the whole DB'):
												
						db_freq = st.slider(
							' ',
							value=(0.0, 1.0),
							max_value=1.0,
							step=0.0001,
							format='%f')
					else: 
						db_freq = None

					if st.checkbox('Filter Unique Individual Count', 
						help='Number of unique individuals in a cluster in the whole DB'):
						count = st.slider(
							'        ',
							value=(1,n_sub),
							max_value=n_sub
							)
					else: 
						count = None

					if st.checkbox('Filter gnomAD Frequency',
						help='gnomAD v2 population frequency (if matched)'):

						# col_gnomad1, col_gnomad2, col_gnomad3 = st.columns(3)

						# with col_gnomad1:
						# 	gnomad_btn_1 = st.button('< 1%')
						# with col_gnomad2:
						# 	gnomad_btn_5 = st.button('< 5%')
						# with col_gnomad3:
						# 	gnomad_btn_10 = st.button('< 10%')
						gnomad_freq_value=(0.0, 1.0)
						# if gnomad_btn_1:
						# 	gnomad_freq_value = (0.0, 0.01)
						# if gnomad_btn_5:
						# 	gnomad_freq_value = (0.0, 0.05)
						# if gnomad_btn_10:
						# 	gnomad_freq_value = (0.0, 0.10)		

						gnomAD_freq = st.slider(
							'',
							value=gnomad_freq_value,
							max_value=1.0,
							step=0.0001,
							format='%f')
					else: 
						gnomAD_freq = None	
					
					if st.checkbox('Filter RefSeq Count',
						help='Number of NCBI curated RefSeq genes overlapping a call'):
						col_rs1, col_rs2 = st.columns(2)
						with col_rs1:
							RefSeq_min = st.number_input(
								'min',
								min_value=0
								)
						with col_rs2:
							RefSeq_max = st.number_input(
								'max',
								min_value=RefSeq_min+1
								)
					else: 
						RefSeq_min = None
						RefSeq_max = None

					if st.checkbox('Filter OMIM Count',
						help='Number of OMIM disease-causing genes overlapping a call'):
						col_omim1, col_omim2 = st.columns(2)
						with col_omim1:
							OMIM_min = st.number_input(
								' min',
								min_value=0
								)
						with col_omim2:
							OMIM_max = st.number_input(
								' max',
								min_value=OMIM_min+1
								)
					else: 
						OMIM_min = None
						OMIM_max = None


					

			num_rows = st.selectbox(
			'Show Number of Rows:',
			[100,500,1000,5000,'All'],
			help='Use lower number of rows to test your query for faster query time')

		qry_dict={
			'chrom1': chrom1,
			'start': start,
			'end': end,
			'svtype': svtype,
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
			'family': family,
			'pt_id': pt_id,
			'num_rows': num_rows,
			'cluster_id': cluster_id,
			'genotype': genotype,
			'disrupt_gene_left': disrupt_gene_left,
			'disrupt_gene_right': disrupt_gene_right,
			'disrupt_repeat_left': disrupt_repeat_left,
			'disrupt_repeat_right': disrupt_repeat_right}

		def build_qry(qry_dict):
			qry = 'SELECT * FROM P2_DB WHERE 1=1'

			if qry_dict['chrom1'] is not None:
				chrom1 = qry_dict['chrom1']
				qry += f' AND chrom1 = "{chrom1}"'

			if qry_dict['start'] is not None:
				start = qry_dict['start']
				qry += f' AND pos1 >= {start}'

			if qry_dict['end'] is not None:
				end = qry_dict['end']
				qry += f' AND pos2 <= {end}'

			if len(qry_dict['svtype']) != 0:
				svtype = '", "'.join(qry_dict['svtype'])
				qry += f' AND SV_type IN ("{svtype}")'

			if qry_dict['sv_len_min'] is not None:
				sv_len_min = qry_dict['sv_len_min']
				qry += f' AND SV_len >= {sv_len_min}'

			if qry_dict['sv_len_max'] is not None:
				sv_len_max = qry_dict['sv_len_max']
				qry += f' AND SV_len <= {sv_len_max}'
			
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

			if len(qry_dict['family']) != 0:
				family = '", "'.join(qry_dict['family'])
				qry += f' AND family IN ("{family}")'

			if len(qry_dict['pt_id']) != 0:
				pt_id = '", "'.join(qry_dict['pt_id'])
				qry += f' AND pt_id IN ("{pt_id}")'

			if qry_dict['genotype'] == 'Homozygous':
				qry += ' AND genotype LIKE "1/1"'
			elif qry_dict['genotype'] == 'Heterozygous':
				qry += ' AND genotype LIKE "0/1"'

			if qry_dict['cluster_id'] is not None:
				cluster_id = qry_dict['cluster_id']
				qry += f' AND cluster = {cluster_id}'

			if qry_dict['disrupt_gene_left']: 
				qry+= f' AND RefSeq_Disrupt_left IS NOT NULL'

			if qry_dict['disrupt_gene_right']: 
				qry+= f' AND RefSeq_Disrupt_right IS NOT NULL'

			if qry_dict['disrupt_repeat_left']: 
				qry+= f' AND RepeatMask_Disrupt_left IS NOT NULL'

			if qry_dict['disrupt_repeat_right']: 
				qry+= f' AND RepeatMask_Disrupt_right IS NOT NULL'

			if qry_dict['num_rows'] == 'All':
				num_rows = 184467440737095516
			else:
				num_rows = qry_dict['num_rows']
			qry += f' LIMIT {num_rows}'
			return qry

		def query():
			qry = build_qry(qry_dict)
			# st.write(db_freq)
			# st.write(qry)
			conn = st.connection('P2_DB', type='sql')
			df=conn.query(qry)
			st.dataframe(
				df,
				hide_index=False,
				height=600)	
			return df

		def plot(df):
			st.title('SV Length Distribution')
			df['Log10_SV_len'] = np.log(df['SV_len'])
			fig = px.histogram(df, x='Log10_SV_len', color='SV_type', marginal = 'box')
			st.plotly_chart(fig, use_container_width=True)
			
			chart1, chart2 = st.columns(2)
			with chart1:
				st.title('SV Type Distribution')
				sv_type_counts = df.SV_type.value_counts()
				fig = px.pie(names=sv_type_counts.index, values=sv_type_counts.values)
				st.plotly_chart(fig)

			with chart2:
				st.title('SD Overlap Distribution')
				SD_counts = df.SD_Overlap.value_counts()
				fig = px.pie(names=SD_counts.index, values=SD_counts.values)
				st.plotly_chart(fig)


		qry_btn = st.button('Query', type='primary')

		nrow_ph = st.empty()


		if qry_btn:
			df=query()
			nrow=df.shape[0]
			nrow_ph.write(f'{nrow} results')
			plot(df)


	with CC_tab2:

		with st.container():

			st.header('Currently querying:')
			st.write('Read depth-based CNV Calls (hg19)')
			col1_cnv, col2_cnv, col3_cnv, col4_cnv, col5_cnv= st.columns(5)

			with col1_cnv:
				chrom1_cnv = st.selectbox(
					'Chromosome ',
					chrom,
					index=None,
					placeholder='Select Chromosome',
					key='chrom_cnv'
					)

				project_cnv = st.multiselect(
				'Project:', 
				projects
				)

				SD_ol_cnv = st.selectbox(
					'SD Overlap ',
					['True', 'False'],
					index=None,
					help = 'If CNV call is overlapping >98% with Segmental Duplications (UCSC)')


			with col2_cnv:
				start_cnv = st.number_input(
					'Left Breakpoint:', 
					value=None, 
					min_value=0)

				is_proband_cnv = st.selectbox(
				'Is CNV call from a proband:',
				['True', 'False'],
				index=None)

				cluster_id_cnv = st.number_input(
					'Cluster ID:', 
					value=None, 
					min_value=-1,
					help='CNV calls are clustered using DBSCAN (eps:500, min_sample:2)')


			with col3_cnv:
				end_cnv = st.number_input(
					'Right breakpoint:', 
					value=None, 
					min_value=0)

				pt_id_cnv = st.multiselect(
				'Individual ID:',
				pt_ids)


				if st.checkbox('Filter CNV Length '):
					col_sv_len1_cnv, col_sv_len2_cnv = st.columns(2)
					with col_sv_len1_cnv:
						sv_len_min_cnv = st.number_input(
							'  min',
							min_value=0
							)
					with col_sv_len2_cnv:
						sv_len_max_cnv = st.number_input(
							'  max',
							min_value=sv_len_min_cnv+1
							)
				else: 
					sv_len_min_cnv = None
					sv_len_max_cnv = None
				

			with col4_cnv:
				svtype_cnv = st.selectbox(
					'CNV Type: ',
					['gain', 'loss'],
					index=None
					)

				family_cnv = st.multiselect(
				'Family ID: ',
				families,
				help='Results will include all family members')




			with col5_cnv:
				level_cnv = st.multiselect(
					'CNV level',
					['DUP', 'TRP', 'MUL_GAIN', 'HET_DEL', 'HOM_DEL', 'UND'],
					help='Characterized level of CNV call')

				OMIM_sym_cnv = st.selectbox(
				'OMIM Disease Genes: ',
				OMIM_syms,
				index=None,
				help='OMIM disease causing genes from UCSC (pheno_key 3/4)'
				)

				rs_sym_cnv = st.selectbox(
				'RefSeq Genes: ',
				rs_syms,
				index=None,
				help='RefSeq genes from NCBI curated list'
				)


			col6_cnv, col7_cnv = st.columns(2)
			with col6_cnv: 
				with st.expander('Breakpoint Filters '):	


					if st.checkbox('Left Breakpoint Disrupting Gene ',
						help='check if left breakpoint is distrupting NCBI curated gene list'):
						disrupt_gene_left_cnv = True
					else: 
						disrupt_gene_left_cnv = False

					if st.checkbox('Right Breakpoint Disrupting Gene ',
						help='check if right breakpoint is distrupting NCBI curated gene list'):
						disrupt_gene_right_cnv = True
					else: 
						disrupt_gene_right_cnv = False


					if st.checkbox('Left Breakpoint In Repeats ',
						help='check if left breakpoint is distrupting RepeatMask repeats'):
						disrupt_repeat_left_cnv = True
					else: 
						disrupt_repeat_left_cnv = False

					if st.checkbox('Right Breakpoint In Repeats ',
						help='check if right breakpoint is distrupting RepeatMask repeats'):
						disrupt_repeat_right_cnv = True
					else: 
						disrupt_repeat_right_cnv = False

					if st.checkbox('CNV overlaps with haploinsufficient genes  ',
						help='HI genes from Collins et al 2022'):
						pHaplo_cnv = True
					else: 
						pHaplo_cnv = False

					if st.checkbox('CNV overlaps with triplosensitive genes  ',
						help='TS genes from Collins et al 2022'):
						pTriplo_cnv = True
					else: 
						pTriplo_cnv = False

			with col7_cnv: 
				with st.expander('Count/Frequency Filters '):	

					if st.checkbox('Filter pseudo-Database Frequency ',
						help='Frequency of number of unique individuals in a cluster in the whole DB'):
												
						db_freq_cnv = st.slider(
							' ',
							value=(0.0, 1.0),
							max_value=1.0,
							step=0.0001,
							format='%f')
					else: 
						db_freq_cnv = None

					if st.checkbox('Filter Unique Individual Count ', 
						help='Number of unique individuals in a cluster in the whole DB'):
						count_cnv = st.slider(
							'        ',
							value=(1,n_sub),
							max_value=n_sub
							)
					else: 
						count_cnv = None

					if st.checkbox('Filter gnomAD Frequency ',
						help='gnomAD v2 population frequency (if matched)'):

						# col_gnomad1, col_gnomad2, col_gnomad3 = st.columns(3)

						# with col_gnomad1:
						# 	gnomad_btn_1 = st.button('< 1%')
						# with col_gnomad2:
						# 	gnomad_btn_5 = st.button('< 5%')
						# with col_gnomad3:
						# 	gnomad_btn_10 = st.button('< 10%')
						gnomad_freq_value=(0.0, 1.0)
						# if gnomad_btn_1:
						# 	gnomad_freq_value = (0.0, 0.01)
						# if gnomad_btn_5:
						# 	gnomad_freq_value = (0.0, 0.05)
						# if gnomad_btn_10:
						# 	gnomad_freq_value = (0.0, 0.10)		

						gnomAD_freq_cnv = st.slider(
							'',
							value=gnomad_freq_value,
							max_value=1.0,
							step=0.0001,
							format='%f')
					else: 
						gnomAD_freq_cnv = None	
					
					if st.checkbox('Filter RefSeq Count ',
						help='Number of NCBI curated RefSeq genes overlapping a call'):
						col_rs1_cnv, col_rs2_cnv = st.columns(2)
						with col_rs1_cnv:
							RefSeq_min_cnv = st.number_input(
								'min',
								min_value=0
								)
						with col_rs2_cnv:
							RefSeq_max_cnv = st.number_input(
								'max',
								min_value=RefSeq_min_cnv+1
								)
					else: 
						RefSeq_min_cnv = None
						RefSeq_max_cnv = None

					if st.checkbox('Filter OMIM Count ',
						help='Number of OMIM disease-causing genes overlapping a call'):
						col_omim1_cnv, col_omim2_cnv = st.columns(2)
						with col_omim1_cnv:
							OMIM_min_cnv = st.number_input(
								' min',
								min_value=0
								)
						with col_omim2_cnv:
							OMIM_max_cnv = st.number_input(
								' max',
								min_value=OMIM_min_cnv+1
								)
					else: 
						OMIM_min_cnv = None
						OMIM_max_cnv = None


					

			num_rows_cnv = st.selectbox(
			'Show Number of Rows: ',
			[100,500,1000,5000,'All'],
			help='Use lower number of rows to test your query for faster query time')

		qry_dict_cnv={
			'chrom1': chrom1_cnv,
			'start': start_cnv,
			'end': end_cnv,
			'svtype': svtype_cnv,
			'level': level_cnv,
			'sv_len_min': sv_len_min_cnv, 
			'sv_len_max': sv_len_max_cnv,
			'SD_ol': SD_ol_cnv, 
			'is_proband': is_proband_cnv,
			'count': count_cnv,
			'db_freq': db_freq_cnv,
			'gnomAD_freq': gnomAD_freq_cnv,
			'project': project_cnv, 
			'OMIM_min': OMIM_min_cnv,
			'OMIM_max': OMIM_max_cnv,
			'RefSeq_min': RefSeq_min_cnv,
			'RefSeq_max': RefSeq_max_cnv,
			'OMIM_sym': OMIM_sym_cnv,
			'rs_sym': rs_sym_cnv,
			'family': family_cnv,
			'pt_id': pt_id_cnv,
			'num_rows': num_rows_cnv,
			'cluster_id': cluster_id_cnv,
			'disrupt_gene_left': disrupt_gene_left_cnv,
			'disrupt_gene_right': disrupt_gene_right_cnv,
			'disrupt_repeat_left': disrupt_repeat_left_cnv,
			'disrupt_repeat_right': disrupt_repeat_right_cnv,
			'pHaplo': pHaplo_cnv,
			'pTriplo': pTriplo_cnv}

		def build_qry_cnv(qry_dict_cnv):
			qry = 'SELECT * FROM CNV_hg19 WHERE 1=1'

			if qry_dict_cnv['chrom1'] is not None:
				chrom1 = qry_dict_cnv['chrom1']
				qry += f' AND chr = "chr{chrom1}"'

			if qry_dict_cnv['start'] is not None:
				start = qry_dict_cnv['start']
				qry += f' AND start >= {start}'

			if qry_dict_cnv['end'] is not None:
				end = qry_dict_cnv['end']
				qry += f' AND end <= {end}'

			if len(qry_dict_cnv['level']) != 0:
				level = '", "'.join(qry_dict_cnv['level'])
				qry += f' AND cnv_lvl IN ("{level}")'

			if qry_dict_cnv['svtype'] is not None:
				svtype = qry_dict_cnv['svtype']
				qry += f' AND type = "{svtype}"'

			if qry_dict_cnv['sv_len_min'] is not None:
				sv_len_min = qry_dict_cnv['sv_len_min']
				qry += f' AND len >= {sv_len_min}'

			if qry_dict_cnv['sv_len_max'] is not None:
				sv_len_max = qry_dict_cnv['sv_len_max']
				qry += f' AND len <= {sv_len_max}'
			
			if qry_dict_cnv['SD_ol'] == 'True':
				SD_ol = 1
				qry += f' AND SD_Overlap = {SD_ol}'
			elif qry_dict_cnv['SD_ol'] == 'False':
				SD_ol = 0
				qry += f' AND SD_Overlap = {SD_ol}'

			if qry_dict_cnv['is_proband'] == 'True':
				is_proband = 1
				qry += f' AND is_proband = {is_proband}'
			elif qry_dict_cnv['is_proband'] == 'False':
				is_proband = 0
				qry += f' AND is_proband = {is_proband}'

			if qry_dict_cnv['count'] is not None:
				count1 = qry_dict_cnv['count'][0]
				count2 = qry_dict_cnv['count'][1]	
				qry += f' AND unique_pt_id_count BETWEEN {count1} and {count2}'

			if qry_dict_cnv['db_freq'] is not None:
				pf1 = qry_dict_cnv['db_freq'][0]
				pf2 = qry_dict_cnv['db_freq'][1]	
				qry += f' AND psuedo_df_freq BETWEEN {pf1} and {pf2}'

			if qry_dict_cnv['gnomAD_freq'] is not None:
				gf1 = qry_dict_cnv['gnomAD_freq'][0]
				gf2 = qry_dict_cnv['gnomAD_freq'][1]	
				qry += f' AND gnomAD_AF BETWEEN {gf1} and {gf2}'

			if len(qry_dict_cnv['project']) != 0:
				project = '", "'.join(qry_dict_cnv['project'])
				qry += f' AND project IN ("{project}")'

			if qry_dict_cnv['OMIM_min'] is not None:
				OMIM_min = qry_dict_cnv['OMIM_min']
				qry += f' AND OMIM_Count >= {OMIM_min}'

			if qry_dict_cnv['OMIM_max'] is not None:
				OMIM_max = qry_dict_cnv['OMIM_max']
				qry += f' AND OMIM_Count <= {OMIM_max}'


			if qry_dict_cnv['RefSeq_min'] is not None:
				RefSeq_min = qry_dict_cnv['RefSeq_min']
				qry += f' AND RefSeq_Count >= {RefSeq_min}'

			if qry_dict_cnv['RefSeq_max'] is not None:
				RefSeq_max = qry_dict_cnv['RefSeq_max']
				qry += f' AND RefSeq_Count <= {RefSeq_max}'


			
			if qry_dict_cnv['OMIM_sym'] is not None:
				OMIM_sym = qry_dict_cnv['OMIM_sym']
				qry += f' AND OMIM_Symbol LIKE "%{OMIM_sym}%"'

			if qry_dict_cnv['rs_sym'] is not None:
				rs_sym = qry_dict_cnv['rs_sym']
				qry += f' AND RefSeq_Symbol LIKE "%{rs_sym}%"'

			if len(qry_dict_cnv['family']) != 0:
				family = '", "'.join(qry_dict_cnv['family'])
				qry += f' AND family IN ("{family}")'

			if len(qry_dict_cnv['pt_id']) != 0:
				pt_id = '", "'.join(qry_dict_cnv['pt_id'])
				qry += f' AND pt_id IN ("{pt_id}")'


			if qry_dict_cnv['cluster_id'] is not None:
				cluster_id = qry_dict_cnv['cluster_id']
				qry += f' AND cluster = {cluster_id}'

			if qry_dict_cnv['disrupt_gene_left']: 
				qry+= f' AND RefSeq_Disrupt_left IS NOT NULL'

			if qry_dict_cnv['disrupt_gene_right']: 
				qry+= f' AND RefSeq_Disrupt_right IS NOT NULL'

			if qry_dict_cnv['disrupt_repeat_left']: 
				qry+= f' AND RepeatMask_Disrupt_left IS NOT NULL'

			if qry_dict_cnv['disrupt_repeat_right']: 
				qry+= f' AND RepeatMask_Disrupt_right IS NOT NULL'

			if qry_dict_cnv['pHaplo']: 
				qry+= f' AND pHaplo_Collins IS NOT NULL'

			if qry_dict_cnv['pTriplo']: 
				qry+= f' AND pTriplo_Collins IS NOT NULL'


			if qry_dict_cnv['num_rows'] == 'All':
				num_rows = 184467440737095516
			else:
				num_rows = qry_dict_cnv['num_rows']
			qry += f' LIMIT {num_rows}'
			return qry

		def query_cnv():
			qry = build_qry_cnv(qry_dict_cnv)
			# st.write(qry)
			conn = st.connection('P2_DB', type='sql')
			df=conn.query(qry)
			st.dataframe(
				df,
				hide_index=False,
				height=600)	
			return df


		qry_btn_cnv = st.button('Query ', type='primary')

		nrow_ph = st.empty()


		if qry_btn_cnv:
			df=query_cnv()
			nrow=df.shape[0]
			nrow_ph.write(f'{nrow} results')


	with CC_tab3:
		with st.container():

			st.header('Currently querying:')
			st.write('Complex Genomic Rearrangement Calls (hg19)')
			col1_cgr, col2_cgr, col3_cgr, col4_cgr, col5_cgr= st.columns(5)

			with col1_cgr:
				chrom1_cgr = st.selectbox(
					'Chromosome ',
					chrom,
					index=None,
					placeholder='Select Chromosome',
					key='chrom_cgr'
					)

				project_cgr = st.multiselect(
				'Project: ', 
				projects
				)

				SD_ol_cgr = st.selectbox(
					'SD Overlap  ',
					['True', 'False'],
					index=None,
					help = 'If CNV call is overlapping >98% with Segmental Duplications (UCSC)')


			with col2_cgr:
				start_cgr = st.number_input(
					'Left Breakpoint: ', 
					value=None, 
					min_value=0)

				is_proband_cgr = st.selectbox(
				'Is CNV call from a proband: ',
				['True', 'False'],
				index=None)

				cluster_id_cgr = st.number_input(
					'Cluster ID: ', 
					value=None, 
					min_value=-1,
					help='CNV calls are clustered using DBSCAN (eps:500, min_sample:2)')


			with col3_cgr:
				end_cgr = st.number_input(
					'Right breakpoint: ', 
					value=None, 
					min_value=0)

				pt_id_cgr = st.multiselect(
				'Individual ID: ',
				pt_ids)


				if st.checkbox('Filter CNV Length  '):
					col_sv_len1_cgr, col_sv_len2_cgr = st.columns(2)
					with col_sv_len1_cgr:
						sv_len_min_cgr = st.number_input(
							'   min',
							min_value=0
							)
					with col_sv_len2_cgr:
						sv_len_max_cgr = st.number_input(
							'   max',
							min_value=sv_len_min_cgr+1
							)
				else: 
					sv_len_min_cgr = None
					sv_len_max_cgr = None
				

			with col4_cgr:
				svtype_cgr = st.selectbox(
					'CNV Type:   ',
					['gain', 'loss'],
					index=None
					)

				family_cgr = st.multiselect(
				'Family ID:  ',
				families,
				help='Results will include all family members')




			with col5_cgr:
				level_cgr = st.multiselect(
					'CNV level ',
					['DUP', 'TRP', 'MUL_GAIN', 'HET_DEL', 'HOM_DEL', 'UND'],
					help='Characterized level of CNV call')

			

			


					

			num_rows_cgr = st.selectbox(
			'Show Number of Rows:  ',
			[100,500,1000,5000,'All'],
			help='Use lower number of rows to test your query for faster query time')

		qry_dict_cgr={
			'chrom1': chrom1_cgr,
			'start': start_cgr,
			'end': end_cgr,
			'svtype': svtype_cgr,
			'level': level_cgr,
			'sv_len_min': sv_len_min_cgr, 
			'sv_len_max': sv_len_max_cgr,
			'SD_ol': SD_ol_cgr, 
			'is_proband': is_proband_cgr,
			'project': project_cgr, 
			'family': family_cgr,
			'pt_id': pt_id_cgr,
			'num_rows': num_rows_cgr,
			'cluster_id': cluster_id_cgr
			}

		def build_qry_cgr(qry_dict_cgr):
			qry = 'SELECT * FROM CGR_hg19 WHERE 1=1'

			if qry_dict_cgr['chrom1'] is not None:
				chrom1 = qry_dict_cgr['chrom1']
				qry += f' AND chr = "chr{chrom1}"'

			if qry_dict_cgr['start'] is not None:
				start = qry_dict_cgr['start']
				qry += f' AND start >= {start}'

			if qry_dict_cgr['end'] is not None:
				end = qry_dict_cgr['end']
				qry += f' AND end <= {end}'

			if len(qry_dict_cgr['level']) != 0:
				level = '", "'.join(qry_dict_cgr['level'])
				qry += f' AND cgr_lvl IN ("{level}")'

			if qry_dict_cgr['svtype'] is not None:
				svtype = qry_dict_cgr['svtype']
				qry += f' AND type = "{svtype}"'

			if qry_dict_cgr['sv_len_min'] is not None:
				sv_len_min = qry_dict_cgr['sv_len_min']
				qry += f' AND len >= {sv_len_min}'

			if qry_dict_cgr['sv_len_max'] is not None:
				sv_len_max = qry_dict_cgr['sv_len_max']
				qry += f' AND len <= {sv_len_max}'
			
			if qry_dict_cgr['SD_ol'] == 'True':
				SD_ol = 1
				qry += f' AND SD_Overlap = {SD_ol}'
			elif qry_dict_cgr['SD_ol'] == 'False':
				SD_ol = 0
				qry += f' AND SD_Overlap = {SD_ol}'

			if qry_dict_cgr['is_proband'] == 'True':
				is_proband = 1
				qry += f' AND is_proband = {is_proband}'
			elif qry_dict_cgr['is_proband'] == 'False':
				is_proband = 0
				qry += f' AND is_proband = {is_proband}'

			if len(qry_dict_cgr['project']) != 0:
				project = '", "'.join(qry_dict_cgr['project'])
				qry += f' AND project IN ("{project}")'

			if len(qry_dict_cgr['family']) != 0:
				family = '", "'.join(qry_dict_cgr['family'])
				qry += f' AND family IN ("{family}")'

			if len(qry_dict_cgr['pt_id']) != 0:
				pt_id = '", "'.join(qry_dict_cgr['pt_id'])
				qry += f' AND pt_id IN ("{pt_id}")'

			if qry_dict_cgr['cluster_id'] is not None:
				cluster_id = qry_dict_cgr['cluster_id']
				qry += f' AND cluster = {cluster_id}'

			qry += ' AND (left_dist <= 1500 OR right_dist <=1500) AND match_CNV_SV != 1'

			if qry_dict_cgr['num_rows'] == 'All':
				num_rows = 184467440737095516
			else:
				num_rows = qry_dict_cgr['num_rows']

			qry += f' LIMIT {num_rows}'
			return qry

		def query_cgr():
			qry = build_qry_cgr(qry_dict_cgr)
			# st.write(qry)
			conn = st.connection('P2_DB', type='sql')
			df=conn.query(qry)
			st.dataframe(
				df,
				hide_index=False,
				height=600)	
			return df


		qry_btn_cgr = st.button('Query  ', type='primary')

		nrow_ph = st.empty()


		if qry_btn_cgr:
			df=query_cgr()
			nrow=df.shape[0]
			nrow_ph.write(f'{nrow} results')