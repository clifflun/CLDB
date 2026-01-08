import pandas as pd
import numpy as np
import streamlit as st
from query.search_filters import search_filters, search_build, igv_filters, vizCNV_filters
from helper.variant_snapshot import variant_snapshot, vizCNV
import streamlit.components.v1 as components

def query(qry, param, db):
	conn = st.connection(db, type='sql')
	df=conn.query(qry, params=param)
	return df   

def final_query(dataset, pipeline, ref, db):
	key=dataset+'_'+pipeline+'_'+ref
	table=pipeline+'_'+ref
	
	if f'{key}_df' not in st.session_state:
		st.session_state[f'{key}_df'] = None

	if f'{key}_displayed' not in st.session_state:
		st.session_state[f'{key}_displayed'] = None

	qry_dict=search_filters(db, key)
	qry, param=search_build(qry_dict, table)
	# st.write(qry) #debug
	if st.button('Query', type='primary', key=f'qry-btn-{key}'):
		df=query(qry, param, db)
		st.session_state[f'{key}_df'] = df	
		st.session_state[f'{key}_displayed']=True

	if st.session_state[f'{key}_displayed']:
		st.markdown("""---""")
		st.write(f'{len(st.session_state[f"{key}_df"].index)} results')
		st.dataframe(st.session_state[f'{key}_df'], hide_index=False, height=600)

	
	if st.session_state[f'{key}_displayed'] and 'CGR' not in key:
		st.markdown("""---""")
		uuids=st.session_state[f'{key}_df'].UUID.unique()
		col1, col2 = st.columns(2)
		with col1:
			uuid=st.selectbox('Select UUID', 
				uuids,
				index=None, 
				key=f'uuid-{key}')
		with col2:
			margin=st.number_input('Margin',
				min_value=0, 
				max_value=10000000,
				value=25000,
				key=f'margin-{key}')
		if st.button('Variant Snapshot', key=f'variant-snapshot-btn-{key}'):
			if uuid is not None:
				qa_df=st.session_state[f'{key}_df'][st.session_state[f'{key}_df']['UUID']==uuid]
				variant_snapshot(qa_df, dataset, pipeline, ref, db, margin)
			else: 
				st.write('Please select UUID')


def igv_query(dataset, pipeline, ref, db):
	key=dataset+'_'+pipeline+'_'+ref
	table=pipeline+'_'+ref

	pt_id = igv_filters(db, key)
	if pt_id:
		if 'SR' in db:
			if db == 'CLDB_SR' and 'hg19' in key: 
				meta = pd.read_csv('/CLDB/meta/meta_SR_hg19.tsv', sep='\t', index_col=False)
			elif db == 'CLDB_SR' and 'hg38' in key: 
				meta = pd.read_csv('/CLDB/meta/meta_SR_hg38.tsv', sep='\t', index_col=False)
			elif db == 'GREGoR_SR': 
				meta = pd.read_csv('/CLDB/meta/gregor_meta.tsv', sep='\t', index_col=False)
			meta = meta[meta['pt_id'] == pt_id]
			bam_path = meta['BAM_path'].squeeze()
		if db == 'CLDB_LR':
			meta = pd.read_csv('/CLDB/meta/meta_LR.tsv', sep='\t', index_col=False)
			meta = meta[meta['pt_id'] == pt_id]
			if ref == 'hg19':
				bam_path = meta['BAM_path_hg19'].squeeze()
			elif ref == 'hg38':
				bam_path = meta['BAM_path_hg38'].squeeze()

		bam_ext = bam_path.split('.')[-1]
		if bam_ext == 'bam':
			bam_index = bam_path+'.bai'
		elif bam_ext == 'cram': 
			bam_index = bam_path+'.crai'

		if ref == 'hg19':
			rmsk_path='/rmsk_hg19_sorted.bed.gz'
		elif ref == 'hg38':
			rmsk_path='/rmsk_hg38_sorted.bed.gz'

		igv_html = f"""
		<!DOCTYPE html>
		<html lang="en">
		<head>
		    <script src="https://cdn.jsdelivr.net/npm/igv@3.2.5/dist/igv.min.js"></script>
		    <style>
		        #igv-container {{
		            width: 100%;
		            height: 600px;
		        }}
		    </style>
		</head>
		<body>
		    <div id="igv-container"></div>
		    <script>
		        document.addEventListener("DOMContentLoaded", function () {{
		            var igvContainer = document.getElementById("igv-container");
		            var options = {{
		                genome: "{ref}",
		                showCenterGuide: true,	                
		                tracks: [
		                    {{
		                        name: "{pt_id}",
		                        url: "http://pnri-app17.pnri.local:3000{bam_path}",
		                        indexURL: "http://pnri-app17.pnri.local:3000{bam_index}",
		                        type: "alignment",
		                        showSoftClips: true
		                    }},
		                    {{
		                    	name: "RepeatMask", 
		                    	type: "annotation", 
		                    	format: "bed", 
		                    	url: "http://pnri-app17.pnri.local:3000{rmsk_path}",
		                    	indexURL: "http://pnri-app17.pnri.local:3000{rmsk_path}.tbi"
		                    }}
		                ]
		            }};
		            igv.createBrowser(igvContainer, options).then(function (browser) {{
		                console.log("IGV.js Ready!"); 
		            }});
		        }});
		    </script>
		</body>
		</html>
		"""
		components.html(igv_html, height=1000)




def vizCNV_query(dataset, pipeline, ref, db):
	key=dataset+'_'+pipeline+'_'+ref
	table=pipeline+'_'+ref

	chrom1, start, end, pt_id = vizCNV_filters(db, key)

	if st.button('Generate vizCNV image (~15-30sec)', key=f'vizCNV-btn-{key}'):
		if db == 'CLDB_SR' and 'hg19' in key: 
			meta = pd.read_csv('/CLDB/meta/meta_SR_hg19.tsv', sep='\t', index_col=False)
		elif db == 'CLDB_SR' and 'hg38' in key: 
			meta = pd.read_csv('/CLDB/meta/meta_SR_hg38.tsv', sep='\t', index_col=False)
		elif db == 'GREGoR_SR': 
			meta = pd.read_csv('/CLDB/meta/gregor_meta.tsv', sep='\t', index_col=False)
		meta = meta[meta['pt_id'] == pt_id]
		P2_path = meta['P2_path'].squeeze()
		MD_path = meta['MD_path'].squeeze()
		GATK_path = meta['GATK_path'].squeeze()
		is_trio = meta['is_trio'].squeeze()
		pr_par2 = P2_path.split(':')[1]
		pr_df = MD_path.split(':')[1]
		pr_seg = pr_df.replace('.regions.bed.gz', '_SLM_segments.tsv')

		if is_trio==1:
			is_trio = "TRUE"
		else: 
			is_trio = "FALSE"

		if 'chr' not in str(chrom1):
			chrom = 'chr' + str(chrom1) 
		else:
			chrom = chrom1
		
		vizCNV(chrom, start, end, ref, pr_df, pr_seg, pr_par2=pr_par2, is_trio=is_trio, gvcf=GATK_path, margin=0, highlight="FALSE")
		st.image(f'/tmp/r_ggplot.png', use_column_width=True)
