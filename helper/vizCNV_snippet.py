
# import rpy2.robjects as robjects
# from rpy2.robjects.packages import importr
# from rpy2.robjects import pandas2ri
# from rpy2.robjects.conversion import localconverter
# import io
# from rpy2.robjects.lib import grdevices


#R.utils, data.table, ggplot2, dplyr, optparse, cowplot, arrow
#if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# VariantAnnotation
#BiocManager::install("GenomicRanges")

#10_127190418_10_127201102_INV0013208SUR_BH10117-1
# BH13731-1 chr12 INV LR hg38

def normalization(rd):
	'''
	Normalize by chr
	'''
	rd['ratio'] = rd.groupby('chrom')['coverage'].transform(lambda x: x / (x.median() + 0.00001))
	return rd

def modify_chrom(value):
	if value.startswith("chr"):
		return value
	else:
		return "chr" + value

def prep_df(df):
	df.columns = ['chrom', 'start', 'end', 'coverage']
	df.iloc[:, 0] = df.iloc[:, 0].astype(str)
	df['chrom'] = df['chrom'].apply(modify_chrom)
	return df

def get_chr_df(input_file, chromosome):
	df = pd.read_csv(input_file, compression='gzip', sep='\t', dtype={0: str})
	df = prep_df(df)
	df = normalization(df)
	df=df[df['chrom'] == chromosome]
	return df

def get_seg_df(input_file, chromosome):
	seg_py=pd.read_csv(input_file, sep="\t", header=None)
	seg_py=seg_py.set_axis(['chrom', 'start', 'end', 'length', 'seg_mean'], axis=1)
	seg_py=seg_py[seg_py['chrom'] == chromosome]
	return seg_py


# pandas2ri.activate()

	# df_file="/Projects/Turkish_ROH/WGS/hg19/processed/mosdepth/BH13674-1_MAPQ30.regions.bed.gz"
	# df_py = get_chr_df(df_file, chrom)
	# df_py = df_py[(df_py['start'] >= from_r-100000) & (df_py['start'] <= to_r+100000)]

	# pr_seg_file='/Projects/Turkish_ROH/WGS/hg19/processed/mosdepth/BH13674-1_MAPQ30_SLM_segments.tsv'
	# mo_seg_file='/Projects/Turkish_ROH/WGS/hg19/processed/mosdepth/BH13674-2_MAPQ30_SLM_segments.tsv'
	# fa_seg_file='/Projects/Turkish_ROH/WGS/hg19/processed/mosdepth/BH13674-3_MAPQ30_SLM_segments.tsv'

	# pr_seg_py=get_seg_df(pr_seg_file, chrom)
	# mo_seg_py=get_seg_df(mo_seg_file, chrom)
	# fa_seg_py=get_seg_df(fa_seg_file, chrom)

	# # Convert Python DataFrame to R DataFrame
	# with localconverter(robjects.default_converter + pandas2ri.converter):
	# 	pr_seg_r = robjects.conversion.py2rpy(pr_seg_py)
	# 	mo_seg_r = robjects.conversion.py2rpy(mo_seg_py)
	# 	fa_seg_r = robjects.conversion.py2rpy(fa_seg_py)
	# 	df_r = robjects.conversion.py2rpy(df_py)
		

	# 	# Import ggplot2
	# 	ggplot2 = importr("ggplot2")

	# 	# Define R code to create a plot using the passed DataFrame
	# 	r_code = """
	# 	function(df, pr_seg, from, to, mo_seg = NULL, fa_seg = NULL) {
	# 		library(ggplot2)
	# 		library(dplyr)

	# 		options(scipen = 999)
	# 		style_rd <- theme_classic()+
	# 	  theme(plot.title = element_text(face = "bold", size = 12),
	# 			legend.position = "top",
	# 			legend.title = element_text(colour="black", size=12),
	# 			legend.text = element_text(size = 12),
	# 			panel.border = element_blank(),
	# 			panel.grid.minor.y = element_blank(),
	# 			panel.grid.minor.x = element_line(linetype = 4,colour = "grey85"),
	# 			panel.grid.major.y = element_line(linetype = 5,colour = "grey70"),
	# 			panel.grid.major.x = element_line(linetype = 5,colour = "grey50"),
	# 			panel.background = element_blank(),
	# 			axis.text.y = element_text(color = "black", size = 30),
	# 			axis.text.x = element_text(color = c("black", "white"), size = 30),
	# 			axis.title = element_text(color = "black", size = 32),
	# 			axis.ticks = element_line(color = "black"))
	# 	scale_rd <- scale_y_continuous(name="Log2 Ratio",
	# 								   limits=c(-2.5, 2),
	# 								   breaks = c(-2,
	# 									 round(log2(1/2),2),
	# 									 round(log2(2/2),2),
	# 									 round(log2(3/2),2),
	# 									 round(log2(4/2),2),
	# 									 round(log2(5/2),2),
	# 									 round(log2(6/2),2)
	# 								   ))

	 
	# 		min.num.mark=100
	# 		bin.size=1000 ## default 1000 bp bin size
	# 		del.log2width <- 0.1
	# 		dup.log2width <- 0.15
	# 		del.range <- c(log2(1*(1-del.log2width )/2),log2(1*(1+del.log2width )/2))
	# 		dup.range <- c(log2(3*(1-dup.log2width )/2),log2(3*(1+dup.log2width )/2))
			
	# 		loss <- pr_seg %>%
	# 			filter(length > min.num.mark) %>% 
	# 			filter(dplyr::between(seg_mean,del.range[1],del.range[2]))

	# 		gain <- pr_seg %>%
	# 			filter(length >= min.num.mark) %>% 
	# 			filter(seg_mean > dup.range[1])


	# 		p <- ggplot(pr_seg,aes(x = end, y = seg_mean))+
	# 		   geom_point(data = subset(df, ratio < 0.7),aes(start,log2(ratio+0.00001)),color="green", size=0.8)+
	#        	   geom_point(data = subset(df, ratio > 1.3),aes(start,log2(ratio+0.00001)),color="red", size=0.8)+
	# 		   geom_segment(data=subset(pr_seg, length>10), aes(x = start, y = seg_mean, xend = end, yend = seg_mean), color="blue", linewidth=1.25)+
	# 		   style_rd+
	# 		   scale_rd+
	# 		   scale_x_continuous(n.breaks=10)+
	# 		   coord_cartesian(expand=F)+
	# 		   labs(x="")+
	# 		   xlim(c(from-50000, to+50000))

	# 		if (!is.null(mo_seg)) {
	# 		p <- p+
	# 		geom_segment(data=subset(mo_seg, length>10), aes(x = start, y = seg_mean, xend = end, yend = seg_mean), color="#E69F00", linewidth=1.25)
	# 		}

	# 		if (!is.null(fa_seg)) {
	# 		p <- p+
	# 		geom_segment(data=subset(fa_seg, length>10), aes(x = start, y = seg_mean, xend = end, yend = seg_mean), color="#39918C", linewidth=1.25)
	# 		}


	# 		# Save the plot
	# 		ggsave("/tmp/r_ggplot.png", plot = p, dpi=300, height=8, width=45, bg="white")
	# 	}
	# 	"""

	# 	# Convert R code into a callable R function
	# 	r_plot_function = robjects.r(r_code)

		# Execute the R function with the Python DataFrame
		# r_plot_function(df_r, pr_seg_r, from_r, to_r, mo_seg_r, fa_seg_r)



	# st.image("/tmp/r_ggplot.png", use_column_width=True)

	# with open("/tmp/r_ggplot.png", "rb") as file:
	# 	st.download_button(
	# 		label="Download image",
	# 		data=file,
	# 		file_name="r_ggplot.png",
	# 		mime="image/png",
	# 	)

