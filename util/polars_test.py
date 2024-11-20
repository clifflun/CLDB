import pandas as pd
import numpy as np
import polars as pl
from time import time 


def timer_func(func): 
    # This function shows the execution time of  
    # the function object passed 
    def wrap_func(*args, **kwargs): 
        t1 = time() 
        result = func(*args, **kwargs) 
        t2 = time() 
        print(f'Function {func.__name__!r} executed in {(t2-t1):.4f}s') 
        return result 
    return wrap_func 




@timer_func
def filter(rm):
	print(rm.filter((pl.col('chrom') == '1') & (pl.col('start') <= 829170) &  (pl.col('end') >= 829170)).to_pandas())

@timer_func
def scan():
	df = (
		pl
		.scan_parquet('../reference/gnomad_v2.1_sv.sites.extracted.parquet')
		# .filter(pl.col('chrom') == '1')
		# .filter(pl.col('start') <= 819170)
		# .filter(pl.col('end') >= 8039170)
		).collect().to_pandas()

	print(df.head(10))

def to_parquet():
	df=pd.read_csv('../reference/gnomad_v4.1_sv.sites.extracted.tsv', sep='\t', dtype={'chrom1':str, 'chrom2':str, 'af':str})
	print(df)
	df.to_parquet('../reference/gnomad_v4.1_sv.sites.extracted.parquet')

def annotate_RepeatMask_pd(chrom1, start, chrom2, end, SV_type, RepeatMask):	

	RepeatMask_splice_left = RepeatMask[(RepeatMask.chrom == chrom1) & (RepeatMask.start<=start) & (RepeatMask.end>=start)]
	RepeatMask_splice_right = RepeatMask[(RepeatMask.chrom == chrom2) & (RepeatMask.start<=end) & (RepeatMask.end>=end)]

	disrupt_left = RepeatMask_splice_left['rep_name'].tolist()
	disrupt_right = RepeatMask_splice_right['rep_name'].tolist()

	return disrupt_left, disrupt_right

def main():
	pd.set_option('display.max_columns', None)
	# to_parquet()
	# scan()
	# print(pl.read_parquet('../reference/gnomad_v4.1_sv.sites.extracted.parquet'))

if __name__ == '__main__':
	main()