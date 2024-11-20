import gzip

infile = 'gnomad.v4.1.sv.sites.vcf.gz'
outfile = 'gnomad_v4.1_sv.sites.extracted.tsv'

with gzip.open(infile, 'rt') as inf, open(outfile, 'wt') as outf:
	for line in inf:
		if line.startswith('#'):
			continue    
		else:
			tmp = line.strip().split('\t')
			chrom1 = tmp[0]
			chrom2 = tmp[0]
			pos1 = tmp[1]
			SV_id = tmp[2]
			for info in tmp[7].split(';'):
				if info.startswith('END='):
					pos2 = info[4:]
				if info.startswith('CHR2'):
					chrom2 = info[5:]
				if info.startswith('SVTYPE='):
					SV_type = info[7:]
				if info.startswith('SVLEN='):
					SV_len = info[6:]
				if info.startswith('AF='):
					af = info[3:]

			out=[chrom1, pos1, chrom2, pos2, SV_type, SV_id, SV_len, af]
			outf.write('\t'.join(out))
			outf.write('\n')