import gzip
import os
import pandas as pd
import numpy as np
import polars as pl
import time
import math

REF_PATH='Z:/Members/clun/CLDB'



def load_SD():
	#get SD data
	row_list=[]
	with open(f'{REF_PATH}/reference/SegDup_hg19_UCSC_no_PAR_merged_1k.bed', 'r') as f: 
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
	OMIM = pd.read_csv(f'{REF_PATH}/reference/OMIM_gene2_hg19_UCSC_all.bed', delimiter='\t', index_col=False)
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
	gnomad=pl.read_parquet(f'{REF_PATH}/reference/gnomad_v2.1_sv.sites.extracted.parquet')
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


def load_decipher():
	decipher=pd.read_csv(f'{REF_PATH}/reference/cnv/population_cnv_grch37_2015_sept.txt', sep='\t')
	print(decipher[decipher['type'] == 0])
	return decipher

def annotate_decipher(chrom, start, end, SV_type, decipher, thresh=500):
	if start == end:
		return [], []

	decipher_splice = decipher.filter(
		(pl.col('svtype') == SV_type)&
		(pl.col('chrom1') == chrom)&
		(abs(pl.col('pos1')-start) < thresh*2)&
		(abs(pl.col('pos2')-end) < thresh*2)
		).to_pandas()

	_id = gnomad_splice['svid'].tolist()
	freq = gnomad_splice['af'].tolist()

	return _id, freq





def load_RefSeq():
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
	rm=pl.read_parquet(f'{REF_PATH}/reference/RepeatMask_hg19.parquet')
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
	IDR = pl.read_parquet(f'{REF_PATH}/reference/cnv/HG19_WG_Collapsed_IDR.parquet')
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
	SD = load_SD()
	OMIM = load_OMIM()
	gnomad = load_gnomAD()
	RS = load_RefSeq()
	RM = load_RepeatMask()
	IDR = load_IDR()
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

		a = annotate_SD(chrom, start, end, SV_type, SD)
		SD_list.append(a)

		sym, dis, inh = annotate_OMIM(chrom, start, end, SV_type, OMIM)	
		OMIM_count.append(len(sym))
		OMIM_symbols.append(';'.join(sym))
		OMIM_disease.append(';'.join(dis))
		OMIM_inh.append(';'.join(inh))

		rs_sym, rs_disrupt_left, rs_disrupt_right = annotate_RefSeq(chrom, start, end, SV_type, RS)	
		RS_count.append(len(rs_sym))
		RS_symbols.append(';'.join(rs_sym))
		RS_disrupt_left.append(';'.join(rs_disrupt_left))
		RS_disrupt_right.append(';'.join(rs_disrupt_right))
		
		rm_disrupt_left, rm_disrupt_right = annotate_RepeatMask(chrom, start, end, SV_type, RM)	
		RM_disrupt_left.append(';'.join(rm_disrupt_left))
		RM_disrupt_right.append(';'.join(rm_disrupt_right))

		idr_disrupt_left, idr_disrupt_right = annotate_IDR(chrom, start, end, SV_type, IDR)	
		IDR_disrupt_left.append(';'.join(idr_disrupt_left))
		IDR_disrupt_right.append(';'.join(idr_disrupt_right))

		gid, af = annotate_gnomad(chrom, start, end, SV_type, gnomad)	
		gnomad_id.append(';'.join(gid))
		gnomad_af.append(';'.join(af))

	annotation={
		'SD_Overlap': SD_list, 
		'OMIM_Count': OMIM_count, 
		'OMIM_Symbol': OMIM_symbols, 
		'OMIM_Disease': OMIM_disease, 
		'OMIM_Inh': OMIM_inh, 
		'RefSeq_Count': RS_count, 
		'RefSeq_Symbol': RS_symbols, 
		'RefSeq_Disrupt_left': RS_disrupt_left,
		'RefSeq_Disrupt_right': RS_disrupt_right, 
		'RepeatMask_Disrupt_left': RM_disrupt_left,
		'RepeatMask_Disrupt_right': RM_disrupt_right, 
		'IDR_Disrupt_left': IDR_disrupt_left,
		'IDR_Disrupt_right': IDR_disrupt_right, 
		'gnomAD_id': gnomad_id, 
		'gnomAD_AF': gnomad_af, 
	}

	df_annotated = df.assign(**annotation)
	df_annotated = df_annotated.drop(columns=['chrom1'])
	return df_annotated


def annotate_CNV_syndromes(chrom, start, end, CNV_syndromes, thresh=500):
	if start == end:
		return []

	CNV_splice = CNV_syndromes[
	    (CNV_syndromes['chrom'] == chrom) &
	    (abs(CNV_syndromes['start'] - start) < thresh * 2) &
	    (abs(CNV_syndromes['end'] - end) < thresh * 2) 
	]

	
	name = CNV_splice['syndrome'].tolist()

	return name

def annotate_step2(df):
	#annotation wrapper
	
	pHaplo = pd.read_csv(f'{REF_PATH}/reference/cnv/pHaplo_filtered_Collins_2022.tsv', sep='\t', index_col=False)
	pHaplo_syms = pHaplo.name
	pTriplo = pd.read_csv(f'{REF_PATH}/reference/cnv/pTriplo_filtered_Collins_2022.tsv', sep='\t', index_col=False)
	pTriplo_syms = pTriplo.name
	CNV_syndromes = pd.read_csv(f'{REF_PATH}/reference/cnv/DECIPHER_CNV_syndrome_hg19.bed', sep='\t', index_col=False)

	pHaplo_Collins=[]
	pTriplo_Collins=[]
	CNV_syndromes_list=[]

	last = ''

	print('Annotation Step2')
	for row in df.itertuples():
		chrom=row.chr
		start=row.start
		end=row.end
		syms = row.RefSeq_Symbol
		if row.Index %1000 ==0:
			print(row.Index, chrom, start, end)
		if not isinstance(syms, float):
			tmp=syms.split(';')
			pHI_out=list(set(tmp).intersection(pHaplo_syms))
			pTS_out=list(set(tmp).intersection(pTriplo_syms))
		else:
			pHI_out=[]
			pTS_out=[]

		pHaplo_Collins.append(';'.join(pHI_out))
		pTriplo_Collins.append(';'.join(pTS_out))
			
		CNV_syndrome = annotate_CNV_syndromes(chrom, start, end, CNV_syndromes)	
		CNV_syndromes_list.append(';'.join(CNV_syndrome))


	annotation={
		
		'pHaplo_Collins': pHaplo_Collins, 
		'pTriplo_Collins': pTriplo_Collins, 
		'DECIPHER_CNV_Syndromes': CNV_syndromes_list
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
	df['chrom1']=df['chr'].apply(remove_chr)
	df=annotate(df)
	df=collapse(df)
	df=annotate_step2(df)
	df.to_csv(output_file, sep='\t', index=False)


if __name__ == '__main__':
	pd.set_option('display.expand_frame_repr', False)
	main('cnv_calls_031324.tsv', 'cnv_annotated.tsv')

