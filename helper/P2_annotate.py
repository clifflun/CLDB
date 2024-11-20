import gzip
import os
import pandas as pd
import numpy as np
import polars as pl
import time

REF_PATH='/store/carvalho/Members/clun/CLDB'

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

def load_DGV():
	DGV=pl.read_parquet(f'{REF_PATH}/reference/dgv_gs.parquet')
	return DGV


def annotate_DGV(chrom, start, end, SV_type, DGV, thresh=500):
	if (SV_type not in ['DUP', 'DEL']) or start == end:
		return [], []


	DGV_splice = DGV.filter(
		(pl.col('DGV_type') == SV_type)&
		(pl.col('chrom') == chrom)&
		(abs(pl.col('start')-start) < thresh*2)&
		(abs(pl.col('end')-end) < thresh*2)&
		(pl.col('start') >= start)&
		(pl.col('end') <= end)
		).to_pandas()

	_id = DGV_splice['DGV_id'].tolist()
	freq = DGV_splice['af'].tolist()

	return _id, freq


def load_RefSeq():
	rs = pd.read_csv(f'{REF_PATH}/reference/RefSeq_hg19.tsv', sep='\t', index_col=False)
	return rs

def annotate_RefSeq(chrom1, start, chrom2, end, SV_type, rs):
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

def annotate_RepeatMask(chrom1, start, chrom2, end, SV_type, RepeatMask):	
	RepeatMask_splice_left = RepeatMask.filter((pl.col('chrom') == chrom1) & (pl.col('start') <= start) &  (pl.col('end') >= start)).to_pandas()
	RepeatMask_splice_right = RepeatMask.filter((pl.col('chrom') == chrom2) & (pl.col('start') <= end) &  (pl.col('end') >= end)).to_pandas()

	disrupt_left = RepeatMask_splice_left['rep_name'].tolist()
	disrupt_right = RepeatMask_splice_right['rep_name'].tolist()

	return disrupt_left, disrupt_right


def annotate(P2):
	#annotation wrapper
	SD = load_SD()
	OMIM = load_OMIM()
	gnomad = load_gnomAD()
	DGV=load_DGV()
	RS = load_RefSeq()
	RM = load_RepeatMask()
	SD_list=[]
	OMIM_count=[]
	OMIM_symbols = []
	OMIM_disease = []
	OMIM_inh = []
	gnomad_id =[]
	gnomad_af =[]
	DGV_id =[]
	DGV_af =[]
	RS_count=[]
	RS_symbols = []
	RS_disrupt_left = []
	RS_disrupt_right = []
	RM_disrupt_left = []
	RM_disrupt_right = []
	last = ''
	for row in P2.itertuples():
		cur = row.pt_id
		if cur != last: 
			last = cur 
			print(last)
		chrom = row.chrom1
		chrom2 = row.chrom2
		start = row.pos1
		end = row.pos2
		SV_type = row.SV_type
		a = annotate_SD(chrom, start, end, SV_type, SD)
		SD_list.append(a)

		sym, dis, inh = annotate_OMIM(chrom, start, end, SV_type, OMIM)	
		OMIM_count.append(len(sym))
		OMIM_symbols.append(';'.join(sym))
		OMIM_disease.append(';'.join(dis))
		OMIM_inh.append(';'.join(inh))

		rs_sym, rs_disrupt_left, rs_disrupt_right = annotate_RefSeq(chrom, start, chrom2, end, SV_type, RS)	
		RS_count.append(len(rs_sym))
		RS_symbols.append(';'.join(rs_sym))
		RS_disrupt_left.append(';'.join(rs_disrupt_left))
		RS_disrupt_right.append(';'.join(rs_disrupt_right))

		rm_disrupt_left, rm_disrupt_right = annotate_RepeatMask(chrom, start, chrom2, end, SV_type, RM)	
		RM_disrupt_left.append(';'.join(rm_disrupt_left))
		RM_disrupt_right.append(';'.join(rm_disrupt_right))

		gid, af = annotate_gnomad(chrom, start, end, SV_type, gnomad)	
		gnomad_id.append(';'.join(gid))
		gnomad_af.append(';'.join(af))

		dgvid, dgvaf = annotate_DGV(chrom, start, end, SV_type, DGV)	
		DGV_id.append(';'.join(dgvid))
		DGV_af.append(';'.join(dgvaf))

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
		'gnomAD_id': gnomad_id, 
		'gnomAD_AF': gnomad_af, 
		'DGV_id': DGV_id, 
		'DGV_AF': DGV_af
	}

	P2_annotated = P2.assign(**annotation)
	return P2_annotated

import argparse

def main(input_file, output_file):
	# Your main functionality here
	print(f"Input file: {input_file}")
	print(f"Output file: {output_file}")

	df = pd.read_csv(input_file, sep='\t', dtype={'chrom1': str, 'chrom2': str})
	
	out = annotate(df)	

	out.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
	# Create the parser
	parser = argparse.ArgumentParser(description="Annotate Parliament2 calls")

	# Add arguments
	parser.add_argument(
		"-i",
		"--input",
		help="Path to the input file.",
		required=True,
		type=str,
	)
	parser.add_argument(
		"-o",
		"--output",
		help="Path to the output file.",
		required=True,
		type=str,
	)

	# Parse the command-line arguments
	args = parser.parse_args()

	# Call the main function with the provided arguments
	main(args.input, args.output)
