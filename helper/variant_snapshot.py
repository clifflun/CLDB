import streamlit as st
import pandas as pd
import gzip
import subprocess


def bamsnap(bam_path, chrom, pos, ref):
	output_file=f'/tmp/{chrom}_{pos}.png'
	command = [
	    "bamsnap",
	    "-bam", bam_path,
	    "-out", output_file,
	    "-pos", f"{chrom}:{pos}",
	    "-margin", "500",
	    "-title", f"{chrom}:{pos}",
	    "-bamplot", "coverage", "read",
	    "-refversion", ref,
	    "-show_soft_clipped",
	    "-read_color_by", "interchrom"
	]
	try:
		subprocess.run(command, check=True, timeout=10)
		return True
	except subprocess.TimeoutExpired:
		st.write('Plotting timeout. Location might contain simple repeats')
		return False

def variant_snapshot(df, dataset, pipeline, ref, db):
	
	chrom1=df['chrom1'].squeeze()
	chrom2=df['chrom2'].squeeze()
	pos1=df['pos1'].squeeze()
	pos2=df['pos2'].squeeze()
	pt_id=df['PT_ID'].squeeze()
	family=df['FAM_ID'].squeeze()
	project=df['PROJECT'].squeeze()
	SV_len=df['SV_LEN'].squeeze()
	SV_type=df['SV_TYPE'].squeeze()

	L_IDR=df['L_IDR'].squeeze()
	R_IDR=df['R_IDR'].squeeze()
	L_repeatmask=df['L_repeatmask'].squeeze()
	R_repeatmask=df['R_repeatmask'].squeeze()
	
	if ref == 'hg38':
		chrom1 = 'chr'+str(chrom1)
		chrom2 = 'chr'+str(chrom2)

	if db == 'CLDB_SR': 
		meta = pd.read_csv('/CLDB/meta/meta_SR_hg19.tsv', sep='\t', index_col=False)
		meta = meta[meta['pt_id'] == pt_id]
		bam_path = meta['BAM_path'].squeeze()
	if db == 'GREGoR_SR': 
		meta = pd.read_csv('/CLDB/meta/gregor_meta.tsv', sep='\t', index_col=False)
		meta = meta[meta['pt_id'] == pt_id]
		bam_path = meta['BAM_path'].squeeze()

	st.write(f'BAM/CRAM source: {bam_path}')
	info1, info2 = st.columns(2)
	with info1:
		st.write(f'Patient: {pt_id}')
		st.write(f'Family: {family}')
		st.write(f'Project: {project}')
		st.write(f'Variant Type: {SV_type}')
		st.write(f'Variant length:{SV_len}')
	with info2:
		st.write(f'L_IDR: {L_IDR}')
		st.write(f'R_IDR: {R_IDR}')
		st.write(f'L_repeatmask: {L_repeatmask}')
		st.write(f'R_repeatmask: {R_repeatmask}')

	st.markdown("""---""")
		
	# Display the image in Streamlit
	col1, col2 = st.columns(2)
	with col1: 
		jct1=bamsnap(bam_path, chrom1, pos1, ref)
		st.write('Junciton 1')
		st.write(f'{chrom1}: {pos1}')
		if jct1:
			st.image(f'/tmp/{chrom1}_{pos1}.png', use_column_width=True)
	with col2:
		jct2=bamsnap(bam_path, chrom2, pos2, ref)
		st.write('Junction 2')
		st.write(f'{chrom2}: {pos2}')
		if jct2:
			st.image(f'/tmp/{chrom2}_{pos2}.png', use_column_width=True)

	