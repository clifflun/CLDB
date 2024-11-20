import gzip
import os
import pandas as pd
import numpy as np
import polars as pl
import time
import math

# REF_PATH='/store/carvalho/Members/clun/CLDB'
REF_PATH='Z:/Members/clun/CLDB'



def load_SD():
	#get SD data
	row_list=[]
	# with open(f'{REF_PATH}/reference/SegDup_hg19_UCSC_no_PAR_merged_1k.bed', 'r') as f: 
	with open(f'{REF_PATH}/reference/SegDup_hg38_UCSC_sorted_merged_1k.bed', 'r') as f: 
		for line in f:
			tmp_dict={}
			tmp = line.strip().split('\t')
			chrom = tmp[0]
			start = int(tmp[1])
			end = int(tmp[2])
			tmp_dict.update({'chrom':chrom, 'start': start, 'end': end})
			row_list.append(tmp_dict)
	df = pd.DataFrame(row_list)
	return df


def annotate_SD(chrom, start, end, SV_type, SD, thresh = 0.98):	
	#True if overlap more than 98% 
	if SV_type == 'BND' or start == end: 
		return False

	#reduce search space
	SD_splice = SD[SD.chrom == chrom]
	segment_len = abs(start-end)
	for row in SD_splice.itertuples():
		#skip if not overlapping
		if row.start>end or row.end<start:
			continue
		min_end = min(end, row.end)
		max_start = max(start, row.start)
		overlap_len = abs(max_start-min_end)
		overlap_per = overlap_len/segment_len
		#98% threshold
		if overlap_per >= thresh:
			return True
	return False


def load_OMIM():
	#get OMIM data
	OMIM = pd.read_csv(f'{REF_PATH}/reference/OMIM_gene2_hg38_UCSC_all.bed', delimiter='\t', index_col=False)
	# OMIM = pd.read_csv(f'{REF_PATH}/reference/OMIM_gene2_hg19_UCSC_all.bed', delimiter='\t', index_col=False)
	OMIM = OMIM[OMIM.pheno_key.notnull()] #only OMIM with known disease
	OMIM['chrom'] = OMIM['chrom'].str.replace('chr', '')
	OMIM.reset_index(inplace=True, drop=True)
	return OMIM

def annotate_OMIM(chrom, start, end, SV_type, OMIM):
	# report anything that overlaps
	if SV_type == 'BND' or start == end: 
		return [],[],[]
	OMIM_splice = OMIM[(OMIM.chrom == chrom) & (OMIM.start<end) & (OMIM.end>start)]

	if OMIM_splice.empty:
		return [],[],[]

	symbol = OMIM_splice['gene_symbol'].tolist()
	name = OMIM_splice['pheno_name'].tolist()
	inh = OMIM_splice['pheno_inh'].astype(str).tolist()
	
	return symbol, name, inh


def load_gnomAD():
	# gnomad=pl.read_parquet(f'{REF_PATH}/reference/gnomad_v2.1_sv.sites.extracted.parquet')
	gnomad=pl.read_parquet(f'{REF_PATH}/reference/gnomad_v4.1_sv.sites.extracted.parquet')
	print(gnomad)
	return gnomad

def annotate_gnomad(chrom, start, end, SV_type, gnomad, thresh=500):
	if (SV_type not in ['DUP', 'DEL', 'INV']) or start == end:
		return [], []

	gnomad_splice = gnomad.filter(
		(pl.col('svtype') == SV_type)&
		(pl.col('chrom1') == chrom)&
		(abs(pl.col('pos1')-start) < thresh*2)&
		(abs(pl.col('pos2')-end) < thresh*2)
		).to_pandas()

	_id = gnomad_splice['svid'].tolist()
	freq = gnomad_splice['af'].tolist()

	return _id, freq


# def load_decipher():
# 	# decipher=pd.read_csv(f'{REF_PATH}/reference/cnv/population_cnv_grch37_2015_sept.txt', sep='\t')
# 	decipher=pd.read_csv(f'{REF_PATH}/reference/cnv/population_cnv_grch37_2015_sept.txt', sep='\t')
# 	print(decipher[decipher['type'] == 0])
# 	return decipher

# def annotate_decipher(chrom, start, end, SV_type, decipher, thresh=500):
# 	if start == end:
# 		return [], []

# 	decipher_splice = decipher.filter(
# 		(pl.col('svtype') == SV_type)&
# 		(pl.col('chrom1') == chrom)&
# 		(abs(pl.col('pos1')-start) < thresh*2)&
# 		(abs(pl.col('pos2')-end) < thresh*2)
# 		).to_pandas()

# 	_id = gnomad_splice['svid'].tolist()
# 	freq = gnomad_splice['af'].tolist()

# 	return _id, freq





def load_RefSeq():
	# rs = pd.read_csv(f'{REF_PATH}/reference/RefSeq_hg38.tsv', sep='\t', index_col=False)
	rs = pd.read_csv(f'{REF_PATH}/reference/RefSeq_hg19.tsv', sep='\t', index_col=False)
	return rs

def annotate_RefSeq(chrom1, start, end, SV_type, rs):
	chrom2=chrom1
	# report anything that overlaps
	if start == end: 
		return [],[],[]
	
	if SV_type != 'BND':
		rs_splice = rs[(rs.chrom == chrom1) & (rs.start<end) & (rs.end>start)]
		if rs_splice.empty:
			return [],[],[]

	rs_splice_left = rs[(rs.chrom == chrom1) & (rs.start<=start) & (rs.end>=start)]
	rs_splice_right = rs[(rs.chrom == chrom2) & (rs.start<=end) & (rs.end>=end)]

	if SV_type != 'BND':
		symbol = rs_splice['gene_id'].tolist()
	else: 
		symbol = []
	disrupt_left = rs_splice_left['gene_id'].tolist()
	disrupt_right = rs_splice_right['gene_id'].tolist()

	return symbol, disrupt_left, disrupt_right


def load_RepeatMask():
	rm=pl.read_parquet(f'{REF_PATH}/reference/RepeatMask_hg38.parquet')
	# rm=pl.read_parquet(f'{REF_PATH}/reference/RepeatMask_hg19.parquet')
	return rm

def annotate_RepeatMask(chrom1, start, end, SV_type, RepeatMask):	
	chrom2 = chrom1
	RepeatMask_splice_left = RepeatMask.filter((pl.col('chrom') == chrom1) & (pl.col('start') <= start) &  (pl.col('end') >= start)).to_pandas()
	RepeatMask_splice_right = RepeatMask.filter((pl.col('chrom') == chrom2) & (pl.col('start') <= end) &  (pl.col('end') >= end)).to_pandas()

	disrupt_left = RepeatMask_splice_left['rep_name'].tolist()
	disrupt_right = RepeatMask_splice_right['rep_name'].tolist()

	return disrupt_left, disrupt_right


def to_parquet(fn, output_file):
	IDR = pd.read_csv(fn, sep='\t', index_col=False)
	IDR.to_parquet(output_file)

def load_IDR():
	IDR = pl.read_parquet(f'{REF_PATH}/reference/cnv/HG38_WG_Collapsed_IDR.parquet')
	# IDR = pl.read_parquet(f'{REF_PATH}/reference/cnv/HG19_WG_Collapsed_IDR.parquet')
	print(IDR)
	return IDR

def annotate_IDR(chrom1, start, end, SV_type, IDR):	
	chrom2 = chrom1
	IDR_splice_left = IDR.filter((pl.col('chrom') == chrom1) & (pl.col('start') <= start) &  (pl.col('end') >= start)).to_pandas()
	IDR_splice_right = IDR.filter((pl.col('chrom') == chrom2) & (pl.col('start') <= end) &  (pl.col('end') >= end)).to_pandas()

	disrupt_left = IDR_splice_left['id'].tolist()
	disrupt_right = IDR_splice_right['id'].tolist()

	return disrupt_left, disrupt_right


def annotate(df):
	#annotation wrapper
	# SD = load_SD()
	# OMIM = load_OMIM()
	# gnomad = load_gnomAD()
	RS = load_RefSeq()
	# RM = load_RepeatMask()
	# IDR = load_IDR()
	SD_list=[]
	OMIM_count=[]
	OMIM_symbols = []
	OMIM_disease = []
	OMIM_inh = []
	gnomad_id =[]
	gnomad_af =[]
	RS_count=[]
	RS_symbols = []
	RS_disrupt_left = []
	RS_disrupt_right = []	
	RM_disrupt_left = []
	RM_disrupt_right = []
	IDR_disrupt_left = []
	IDR_disrupt_right = []
	last = ''

	print('Annotation Step1')
	for row in df.itertuples():
		# cur = row.pt_id
		# if cur != last: 
		# 	last = cur 
		# 	print(last)
		chrom = row.chrom1
		start = row.start
		end = row.end
		if row.type == 'gain':
			SV_type = 'DUP'
		elif row.type == 'loss':
			SV_type = 'DEL'

		if row.Index %1000 ==0:
			print(row.Index, chrom, start, end, SV_type)

		# a = annotate_SD(chrom, start, end, SV_type, SD)
		# SD_list.append(a)

		# sym, dis, inh = annotate_OMIM(chrom, start, end, SV_type, OMIM)	
		# OMIM_count.append(len(sym))
		# OMIM_symbols.append(';'.join(sym))
		# OMIM_disease.append(';'.join(dis))
		# OMIM_inh.append(';'.join(inh))

		rs_sym, rs_disrupt_left, rs_disrupt_right = annotate_RefSeq(chrom, start, end, SV_type, RS)	
		RS_count.append(len(rs_sym))
		RS_symbols.append(';'.join(rs_sym))
		RS_disrupt_left.append(';'.join(rs_disrupt_left))
		RS_disrupt_right.append(';'.join(rs_disrupt_right))
		
		# rm_disrupt_left, rm_disrupt_right = annotate_RepeatMask(chrom, start, end, SV_type, RM)	
		# RM_disrupt_left.append(';'.join(rm_disrupt_left))
		# RM_disrupt_right.append(';'.join(rm_disrupt_right))

		# idr_disrupt_left, idr_disrupt_right = annotate_IDR(chrom, start, end, SV_type, IDR)	
		# IDR_disrupt_left.append(';'.join(idr_disrupt_left))
		# IDR_disrupt_right.append(';'.join(idr_disrupt_right))

		# gid, af = annotate_gnomad(chrom, start, end, SV_type, gnomad)	
		# gnomad_id.append(';'.join(gid))
		# gnomad_af.append(';'.join(af))

	annotation={
		# 'SD_Overlap': SD_list, 
		# 'OMIM_Count': OMIM_count, 
		# 'OMIM_Symbol': OMIM_symbols, 
		# 'OMIM_Disease': OMIM_disease, 
		# 'OMIM_Inh': OMIM_inh, 
		'RefSeq_Count': RS_count, 
		'RefSeq_Symbol': RS_symbols, 
		'RefSeq_Disrupt_left': RS_disrupt_left,
		'RefSeq_Disrupt_right': RS_disrupt_right, 
		# 'RepeatMask_Disrupt_left': RM_disrupt_left,
		# 'RepeatMask_Disrupt_right': RM_disrupt_right, 
		# 'IDR_Disrupt_left': IDR_disrupt_left,
		# 'IDR_Disrupt_right': IDR_disrupt_right, 
		# 'gnomAD_id': gnomad_id, 
		# 'gnomAD_AF': gnomad_af, 
	}

	df_annotated = df.assign(**annotation)
	df_annotated = df_annotated.drop(columns=['chrom1'])
	return df_annotated


def annotate_overlap(chrom, start, end, df):
	if start == end:
		return []

	#any overlap
	df_splice = df[
	    (df['chrom'] == chrom) &
	    (df['end'] > start) & (df['start'] < end) 
	]
	#get percentage overlap
	pol_seg=[]
	pol_cnv=[]
	if len(df_splice) != 0:
		for row in df_splice.itertuples():
			# print(row.start, row.end)
			# print(chrom, start, end)
			segment_len=row.end-row.start
			min_end = min(end, row.end)
			max_start = max(start, row.start)
			overlap_len = abs(max_start-min_end)

			## %age overlap for db seg
			overlap_per_seg = overlap_len/segment_len
			if overlap_per_seg >=1:
				overlap_per_seg = 1
			pol_seg.append(str(overlap_per_seg))

			## %age overlap for input cnv
			overlap_per_cnv = overlap_len/(end-start)
			if overlap_per_cnv >=1:
				overlap_per_cnv = 1
			pol_cnv.append(str(overlap_per_cnv))
	name = df_splice['syndrome'].tolist()
	# print(name, pol_seg, pol_cnv)
	return name, pol_seg, pol_cnv

def annotate_clinGen(chrom, start, end, df):
	if start == end:
		return []

	# completely encompass
	df_splice = df[
	    (df['chrom'] == chrom) &
	    (df['start'] > start) &
	    (df['end'] < end) 
	]
	
	name = df_splice['symbol'].tolist()

	return name

def annotate_step2(df):
	#annotation wrapper
	
	# pHaplo = pd.read_csv(f'{REF_PATH}/reference/cnv/pHaplo_filtered_Collins_2022.tsv', sep='\t', index_col=False)
	# pHaplo_syms = pHaplo.name
	# pTriplo = pd.read_csv(f'{REF_PATH}/reference/cnv/pTriplo_filtered_Collins_2022.tsv', sep='\t', index_col=False)
	# pTriplo_syms = pTriplo.name
	# CNV_syndromes = pd.read_csv(f'{REF_PATH}/reference/cnv/DECIPHER_CNV_syndrome_hg19.bed', sep='\t', index_col=False)
	# ISCAs = pd.read_csv(f'{REF_PATH}/reference/cnv/ISCA_filtered_hg19.bed', sep= '\t', index_col=False)
	# clinGenHaplo = pd.read_csv(f'{REF_PATH}/reference/cnv/clinGenHaplo_filtered.bed', sep='\t', index_col=False)
	# clinGenHaplo_syms=clinGenHaplo.symbol
	# clinGenTriplo = pd.read_csv(f'{REF_PATH}/reference/cnv/clinGenTriplo_filtered.bed', sep='\t', index_col=False)
	# clinGenTriplo_syms=clinGenTriplo.symbol
	# pHaplo = pd.read_csv(f'{REF_PATH}/reference/cnv/pHaplo_filtered_Collins_2022_hg38.tsv', sep='\t', index_col=False)
	# pHaplo_syms = pHaplo.name
	# pTriplo = pd.read_csv(f'{REF_PATH}/reference/cnv/pTriplo_filtered_Collins_2022_hg38.tsv', sep='\t', index_col=False)
	# pTriplo_syms = pTriplo.name
	CNV_syndromes = pd.read_csv(f'{REF_PATH}/reference/cnv/DECIPHER_CNV_syndrome_hg38.bed', sep='\t', index_col=False)
	ISCAs = pd.read_csv(f'{REF_PATH}/reference/cnv/ISCA_filtered_hg38.bed', sep= '\t', index_col=False)
	# clinGenHaplo = pd.read_csv(f'{REF_PATH}/reference/cnv/clinGenHaplo_filtered_hg38.bed', sep='\t', index_col=False)
	# clinGenTriplo = pd.read_csv(f'{REF_PATH}/reference/cnv/clinGenTriplo_filtered_hg38.bed', sep='\t', index_col=False)


	pHaplo_Collins=[]
	pTriplo_Collins=[]
	CNV_syndromes_list=[]
	CNV_syndromes_pol_seg=[]
	CNV_syndromes_pol_cnv=[]
	ISCA_list=[]
	ISCA_pol_seg=[]
	ISCA_pol_cnv=[]
	clinGenHaplo_list=[]
	clinGenTriplo_list=[]

	last = ''

	print('Annotation Step2')
	for row in df.itertuples():
		chrom=row.chr
		start=row.start
		end=row.end
		syms = row.RefSeq_Symbol
		if row.Index %1000 ==0:
			print(row.Index, chrom, start, end)
		# if not isinstance(syms, float):
		# 	tmp=syms.split(';')
		# 	pHI_out=list(set(tmp).intersection(pHaplo_syms))
		# 	pTS_out=list(set(tmp).intersection(pTriplo_syms))
		# 	clinGen_HI_out=list(set(tmp).intersection(clinGenHaplo_syms))
		# 	clinGen_pTS_out=list(set(tmp).intersection(clinGenTriplo_syms))
		# else:
		# 	pHI_out=[]
		# 	pTS_out=[]
		# 	clinGen_HI_out=[]
		# 	clinGen_pTS_out=[]

		# pHaplo_Collins.append(';'.join(pHI_out))
		# pTriplo_Collins.append(';'.join(pTS_out))
		# clinGenHaplo_list.append(';'.join(clinGen_HI_out))
		# clinGenTriplo_list.append(';'.join(clinGen_pTS_out))

		CNV_syndrome, pol_seg, pol_cnv = annotate_overlap(chrom, start, end, CNV_syndromes)	
		CNV_syndromes_list.append(';'.join(CNV_syndrome))
		CNV_syndromes_pol_seg.append(';'.join(pol_seg))
		CNV_syndromes_pol_cnv.append(';'.join(pol_cnv))


		ISCA, pol_seg, pol_cnv = annotate_overlap(chrom, start, end, ISCAs)	
		ISCA_list.append(';'.join(ISCA))
		ISCA_pol_seg.append(';'.join(pol_seg))
		ISCA_pol_cnv.append(';'.join(pol_cnv))

		

	annotation={
		
		# 'pHaplo_Collins': pHaplo_Collins, 
		# 'pTriplo_Collins': pTriplo_Collins, 
		# 'pHaplo_clinGen': clinGenHaplo_list,
		# 'pTriplo_clinGen': clinGenTriplo_list
		'DECIPHER_CNV_Syndromes': CNV_syndromes_list,
		'DECIPHER_CNV_pol_seg': CNV_syndromes_pol_seg,
		'DECIPHER_CNV_pol_cnv': CNV_syndromes_pol_cnv,
		'ISCA': ISCA_list,
		'ISCA_pol_seg': ISCA_pol_seg,
		'ISCA_pol_cnv': ISCA_pol_cnv
	}

	df_annotated = df.assign(**annotation)
	return df_annotated


def remove_chr(value):
    if value.startswith("chr"):
        return value.split('chr')[1]
    else:
        return value


def collapse(df):
	print('Collapsing calls in repeats')
	grouped = df.groupby(['type', 'IDR_Disrupt_left', 'IDR_Disrupt_right'])['cluster'].agg(['unique'])
	grouped = grouped[grouped['unique'].apply(len) > 1]
	for row in grouped.itertuples():
		#avoid grouping calls without IDR disrupt
		if row[0][1] == '' and row[0][1] == '':
			continue
		print(row)
		idx_ = df[(df['type'] == row[0][0]) & (df['IDR_Disrupt_left'] == row[0][1]) & (df['IDR_Disrupt_right'] == row[0][2])].index
		new_id = df.loc[idx_]['cluster'].max()
		df.loc[idx_, 'cluster'] = new_id


	##recalculate freq
	#get count and propensity
	print('Getting cluster count')
	total_pt = len(df['pt_id'].unique())
	df['count'] = df.groupby(['cluster'])['cluster'].transform('count')
	df.loc[df['cluster'] == -1, 'count'] = 1
	df['cluster_propensity'] = df['count']/total_pt
	
	#get pseudo db freq
	print('Getting pseudo db freq')
	df['unique_pt_id_count'] = df.groupby('cluster')['pt_id'].transform('nunique')
	df['unique_pt_id_count'] = np.where(df['cluster'] == -1, 1, df['unique_pt_id_count'])
	nsub=df.pt_id.nunique()
	df['psuedo_df_freq'] = df['unique_pt_id_count']/nsub
	return df



def main(input_file, output_file):
	df = pd.read_csv(input_file, sep='\t')
	# df['chrom1']=df['chr'].apply(remove_chr)
	# df=annotate(df)	
	# df=collapse(df)
	df=annotate_step2(df)
	df=df[['chr','DECIPHER_CNV_Syndromes','DECIPHER_CNV_pol_seg', 'DECIPHER_CNV_pol_cnv','ISCA','ISCA_pol_seg', 'ISCA_pol_cnv']]
	# df=df[['chr','RefSeq_Count', 'RefSeq_Symbol', 'RefSeq_Disrupt_left','RefSeq_Disrupt_right']]
	df.to_csv(output_file, sep='\t', index=False)
	# CNV_syndromes = pd.read_csv(f'{REF_PATH}/reference/cnv/DECIPHER_CNV_syndrome_hg38.bed', sep='\t', index_col=False)
	# annotate_overlap('chrX', 51996000, 52025000, CNV_syndromes)


import argparse



if __name__ == "__main__":
	# pd.set_option('display.expand_frame_repr', False)
	# CNV_syndromes = pd.read_csv(f'{REF_PATH}/reference/cnv/DECIPHER_CNV_syndrome_hg19.bed', sep='\t', index_col=False)
	# print(CNV_syndromes.head(10))
	# CNV_syndrome = annotate_CNV_syndromes("chr15", 25_165_000, 25_438_000, CNV_syndromes)
	# print(CNV_syndrome)
	main('gregor_cnv_calls_no_outliers_collapsed_070824.tsv', 'test.tsv')
	# g=load_gnomAD()
	# a,b=annotate_gnomad('1', 10000, 295666, 'DUP', g)
	# print(a,b)
	# Create the parser
	# parser = argparse.ArgumentParser(description="Annotate Parliament2 calls")

	# # Add arguments
	# parser.add_argument(
	# 	"-i",
	# 	"--input",
	# 	help="Path to the input file.",
	# 	required=True,
	# 	type=str,
	# )
	# parser.add_argument(
	# 	"-o",
	# 	"--output",
	# 	help="Path to the output file.",
	# 	required=True,
	# 	type=str,
	# )

	# # Parse the command-line arguments
	# args = parser.parse_args()

	# # Call the main function with the provided arguments
	# main(args.input, args.output)
