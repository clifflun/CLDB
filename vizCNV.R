library(ggplot2)
library(dplyr)
library(data.table)
library(optparse)
library(cowplot)

option_list <- list(
  make_option(c("--chrom"), type = "character", default = NA, help = "chrom:from-to 1"),
  make_option(c("--from"), type = "integer", default = NA, help = "chrom:from-to 2"),
  make_option(c("--to"), type = "integer", default = NA, help = "chrom:from-to 3"),
  make_option(c("--margin"), type = "integer", default = NA, help = "window margin"),
  make_option(c("--highlight"), type = "logical", default = NA, help = "highlight"),
  make_option(c("--refver"), type = "character", default = NA, help = "hg19/hg38"),  
  make_option(c("--pr_par2"), type = "character", default = NA, help = "Path to proband parliamen2 input file"),
  make_option(c("--pr_df"), type = "character", default = NA, help = "Path to proband data input file"),
  make_option(c("--pr_seg"), type = "character", default = NA, help = "Path to proband segment input file"),
  make_option(c("--mo_seg"), type = "character", default = NA, help = "Path to mother segment input file"),
  make_option(c("--fa_seg"), type = "character", default = NA, help = "Path to father segment input file"),
  make_option(c("--is_trio"), type = "logical", default = FALSE, help = "trio or not"),
  make_option(c("--gvcf"), type = "character", default = NA, help = "Path to GATK input file")

)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
print(opt)


get_chr_segment <- function(path, chr){
	df <- data.table::fread(path)
	names(df) <- c("chrom", "start", "end", "length", "seg_mean")
	df <- df %>%
		filter(chrom == chr)
	return(df)
}

normalization = function(rd){
  out <- rd %>%
  group_by(chrom) %>%
    dplyr::mutate(ratio=coverage/median(coverage+0.00001))
  return(out)
}

chr <- opt$chrom
from <- opt$from
to <- opt$to
margin <- opt$margin
refver <- opt$refver
pr_par2_path <- opt$pr_par2
pr_df_path <- opt$pr_df
pr_seg_path <- opt$pr_seg
mo_seg_path <- opt$mo_seg
fa_seg_path <- opt$fa_seg
is_trio <- opt$is_trio
highlight <- opt$highlight
gvcf <- opt$gvcf
if (grepl("vcf.gz", opt$gvcf)){
	gvcf <- opt$gvcf
} else {
	gvcf <- NA
}

if (refver == "hg19"){
	ref_path = "/CLDB/util/BEDanno/reference/hg19/"
	} else if (refver == "hg38"){
	ref_path = "/CLDB/util/BEDanno/reference/hg38/"
	}

#get df/seg data
pr_df <- data.table::fread(pr_df_path)
names(pr_df) <- c("chrom", "start", "end", "coverage")
pr_df <- pr_df %>%
			mutate(chrom = ifelse(grepl("^chr", chrom), chrom, paste0("chr", chrom)))
pr_df <- normalization(pr_df)
pr_df <- pr_df %>%
			filter(chrom == chr)
pr_seg <- get_chr_segment(pr_seg_path, chr)

if (!is.na(mo_seg_path)){
	mo_seg <- get_chr_segment(mo_seg_path, chr)
}
if (!is.na(fa_seg_path)){
	fa_seg <- get_chr_segment(fa_seg_path, chr)
}

# boundaries
final_from <- from-margin
if (final_from < 0){
	final_from <- 0
}
final_to <- to+margin
if (final_to > max(pr_df$end)){
	final_to <- max(pr_df$end)
}


style_anno <- theme_classic()+
  theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.title.x = element_blank(),
        panel.grid.major.x = element_line(linetype = 5,colour = "grey50"),
        axis.text.y=element_text(color = "white", size =30),  #remove y axis labels
        axis.ticks.y=element_blank(), 
        axis.title.y = element_text(color = "black", size = 32)
  )

scale_anno <- scale_y_continuous(limits = c(-0.01,.11))


process_sv <- function(sv){
  sv <- sv %>% 
    filter(AVGLEN > 10000 & AVGLEN < 100000000,
           SVTYPE != "BND")
  sv <- sv %>% 
    mutate(colors = case_when(SVTYPE == "DEL" ~ "darkblue",
                             SVTYPE == "DUP" ~ "#8b0000",
                             SVTYPE == "INS" ~ "darkgreen", 
                             SVTYPE == "INV" ~ "darkorange", 
                             TRUE ~ "white")) %>% 
    mutate(idx = sample.int(size = dim(sv)[1], n = 980, replace = T)/10000)
  sv <- sv %>% 
    mutate(start = POS, 
           end = as.integer(END)) %>% 
    relocate(CHROM, start, end) 
  return(sv)
}

chr_nochr <- gsub("chr", "", chr)

p2 <- bedr::read.vcf(pr_par2_path, split.info = T)$vcf
pr_sv <- process_sv(p2)
pr_sv <- pr_sv %>% 
  filter(CHROM == chr)
p_par2 <- ggplot(pr_sv, aes(x = start, y = idx)) +
  annotate("rect", xmin = pr_sv$start, xmax = pr_sv$end, ymin = pr_sv$idx, ymax = pr_sv$idx+0.0001, color = pr_sv$color)+
  style_anno+
  scale_anno+
  ylab("Parliament2")+
	scale_x_continuous(n.breaks=10)+
	coord_cartesian(xlim=c(final_from, final_to), expand=F)


segdup <- fread(paste0(ref_path, "ucsc_segdup.tsv"))
segdup <- segdup %>%
	filter(chrom == chr)
segdup <- segdup %>%
mutate(idx = sample(1:100, size = dim(segdup)[1], replace = T)/1000) %>%
mutate(color = case_when(level == "SD_low" ~ "gray40",
                         level == "SD_mid" ~ "yellow2",
                         level == "SD_high" ~ "darkorange"))
segdup_pos <- segdup %>%
filter(strand == "+")
segdup_neg <- segdup %>%
filter(strand == "-")

p_segdup <- ggplot(segdup, aes(x = start, y = idx)) +
geom_segment(data = segdup_pos, aes(x = start, y = idx, xend = end, yend = idx), arrow = arrow(length = unit(0.05, "inches")), color = segdup_pos$color)+
geom_segment(data = segdup_neg, aes(x = end, y = idx, xend = start, yend = idx), arrow = arrow(length = unit(0.05, "inches")), color = segdup_neg$color)+
style_anno+
scale_anno+
ylab("SegDup")+
scale_x_continuous(n.breaks=10)+
coord_cartesian(xlim=c(final_from, final_to), expand=F)

OMIM <- data.table::fread(paste0(ref_path, "OMIM_sorted.bed"))
names(OMIM) <- c("chrom", "start", "end", "symbol", "disease", "inheritance")
OMIM <- OMIM %>%
filter(chrom == chr_nochr) 
OMIM <- OMIM %>%
mutate(idx = sample(1:100, size = dim(OMIM)[1], replace = T)/1000) %>%
mutate(color = "green4")

p_OMIM <- ggplot(OMIM, aes(x = start, y = idx)) +
annotate("rect", xmin = OMIM$start, xmax = OMIM$end, ymin = OMIM$idx, ymax = OMIM$idx+0.0001, color = OMIM$color)+
annotate("text", x = OMIM$start, y = OMIM$idx, label = OMIM$symbol, size = 6.5, hjust = 1.1)+
style_anno+
scale_anno+
ylab("OMIM")+
scale_x_continuous(n.breaks=10)+
coord_cartesian(xlim=c(final_from, final_to), expand=F)

options(scipen = 999)
style_rd <- theme_classic()+
theme(plot.title = element_text(face = "bold", size = 12),
								legend.position = "top",
								legend.title = element_text(colour="black", size=12),
								legend.text = element_text(size = 12),
								panel.border = element_blank(),
								panel.grid.minor.y = element_blank(),
								panel.grid.minor.x = element_line(linetype = 4,colour = "grey85"),
								panel.grid.major.y = element_line(linetype = 5,colour = "grey70"),
								panel.grid.major.x = element_line(linetype = 5,colour = "grey50"),
								panel.background = element_blank(),
								axis.text.y = element_text(color = "black", size = 30),
								axis.text.x = element_text(color = c("black", "white"), size = 30),
								axis.title = element_text(color = "black", size = 32),
								axis.ticks = element_line(color = "black"))
scale_rd <- scale_y_continuous(name="Log2 Ratio",
					   limits=c(-2.5, 2),
					   breaks = c(-2,
						 round(log2(1/2),2),
						 round(log2(2/2),2),
						 round(log2(3/2),2),
						 round(log2(4/2),2),
						 round(log2(5/2),2),
						 round(log2(6/2),2)
	))


min.num.mark=100
bin.size=1000 ## default 1000 bp bin size
del.log2width <- 0.1
dup.log2width <- 0.15
del.range <- c(log2(1*(1-del.log2width )/2),log2(1*(1+del.log2width )/2))
dup.range <- c(log2(3*(1-dup.log2width )/2),log2(3*(1+dup.log2width )/2))

loss <- pr_seg %>%
filter(length > min.num.mark) %>% 
filter(dplyr::between(seg_mean,del.range[1],del.range[2]))

gain <- pr_seg %>%
filter(length >= min.num.mark) %>% 
filter(seg_mean > dup.range[1])

p0 <- ggplot(subset(pr_df, ratio>=0.7 & ratio <= 1.3),aes(x = end, y = seg_mean))+
		geom_point(aes(start, log2(ratio+0.00001)), color="black")+
		geom_point(data = subset(pr_df, ratio < 0.7),aes(start,log2(ratio+0.00001)),color="green", size=1)+
		geom_point(data = subset(pr_df, ratio > 1.3),aes(start,log2(ratio+0.00001)),color="red", size=1)+
		geom_segment(data=subset(pr_seg, length>10), aes(x = start, y = seg_mean, xend = end, yend = seg_mean), color="blue", linewidth=1.25)+
		style_rd+
		scale_rd+
		scale_x_continuous(n.breaks=10)+
		coord_cartesian(expand=F)+
		labs(x="")

if (isTRUE(highlight)){
	p0 <- p0 +
	annotate("rect", fill = "blue", alpha=0.3, xmin = from, xmax = to, ymin = -Inf, ymax = Inf)		
}

p1 <- ggplot()+
		geom_segment(aes(x=final_from, y=0, xend=0, yend=-1), color="red", linetype="dashed", size=2)+
		geom_segment(aes(x=final_to, y=0, xend=max(pr_df$end), yend=-1), color="red", linetype="dashed", size=2)+
		theme_classic()+
		theme(
			axis.line = element_blank(),
	        axis.text.x=element_blank(),
	        axis.ticks.x=element_blank(), #remove x axis ticks
	        axis.title.x = element_blank(),
	        panel.grid.major.x = element_line(linetype = 5,colour = "white"),
	        axis.text.y=element_text(color = "white", size =30),  #remove y axis labels
	        axis.ticks.y=element_blank(), 
	        axis.title.y = element_text(color = "white", size = 32)
		)+
		coord_cartesian(expand=F)+
		xlim(0,max(pr_df$end))+
		ylim(-1,0)

p <- ggplot(pr_df,aes(x = end, y = seg_mean))+
		geom_point(aes(start, log2(ratio+0.00001)), color="black")+
		geom_point(data = subset(pr_df, ratio < 0.7),aes(start,log2(ratio+0.00001)),color="green", size=1)+
		geom_point(data = subset(pr_df, ratio > 1.3),aes(start,log2(ratio+0.00001)),color="red", size=1)+
		geom_segment(data=subset(pr_seg, length>10), aes(x = start, y = seg_mean, xend = end, yend = seg_mean), color="blue", linewidth=1.25)+
		style_rd+
		scale_rd+
		scale_x_continuous(n.breaks=10)+
		coord_cartesian(xlim=c(final_from, final_to), expand=F)+
		labs(x="")

if (isTRUE(highlight)){
	p <- p +
	annotate("rect", fill = "blue", alpha=0.3, xmin = from, xmax = to, ymin = -Inf, ymax = Inf)		
}

if (!is.na(mo_seg_path)) {
p <- p+
geom_segment(data=subset(mo_seg, length>10), aes(x = start, y = seg_mean, xend = end, yend = seg_mean), color="#E69F00", linewidth=1.25)
}

if (!is.na(fa_seg_path)) {
p <- p+
geom_segment(data=subset(fa_seg, length>10), aes(x = start, y = seg_mean, xend = end, yend = seg_mean), color="#39918C", linewidth=1.25)
}

### baf
style_snp <- theme_classic()+
  theme(
  		# plot.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        legend.title = element_text(colour="black", size=32),
        legend.text = element_text(size = 32),
        panel.border = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_line(linetype = 4,colour = "grey85"),
        panel.grid.major.y = element_line(linetype = 5,colour = "grey70"),
        panel.grid.major.x = element_line(linetype = 5,colour = "grey50"),
        panel.background = element_blank(),
        # axis.line = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.title.x = element_blank(),
        axis.text.y=element_text(color = "black", size =27),  #remove y axis labels
        axis.ticks.y=element_blank(), 
        axis.title.y = element_text(color = "black", size = 32, margin = margin(r = 10))
        )


scale_snp <- scale_y_continuous(name="B-allele frequency",
                                limits = c(-0.05, 1.05),
                                breaks = c(round(0,3),
                                           round(1/2,3),
                                           round(1/3,3),
                                           round(2/3,3),
                                           round(1/4,3),
                                           round(3/4,3),
                                           # round(2/5,3),
                                           # round(3/5,3),
                                           1))


ReadVCF <- function(file, gr, refver){
	vcf <-  VariantAnnotation::readVcf(file = file, genome = refver, param = gr)
	vcf.gr <- vcf@rowRanges
	GT <- VariantAnnotation::geno(vcf)$GT
	AD <- VariantAnnotation::geno(vcf)$AD
	DP <- VariantAnnotation::geno(vcf)$DP
	PR_ID=colnames(GT)[1]
	P1_ID=colnames(GT)[2]
	P2_ID=colnames(GT)[3]
	G1=c('0/0',"0|0")
	G2=c('1/1',"1|1")
	G3=c('0/1',"0|1")
	GT <- as.data.table(GT)
	setnames(GT,colnames(GT),c("index","P1","P2"))
	GT.anno <- GT %>% 
	mutate(B_InhFrom=case_when(index %in% G3 & P1 %in% G1 & P2  %in% c(G2,G3) ~ P2_ID, #cases 11,12
	                           index %in% G3 & P1 %in% c(G2,G3) & P2  %in% G1 ~ P1_ID, #case 13,16
	                           index %in% G2 & P1 %in% G1 & P2  %in% c(G2,G3) ~ P2_ID, #cases 20,21
	                           index %in% G2 & P1 %in% c(G2,G3) & P2  %in% G1 ~ P1_ID, #case 22,25
	                           TRUE ~ "Notphased")) %>% 
	mutate(A_InhFrom=case_when(index %in% G1 & P1 %in% c(G1,G3) & P2  %in% G2 ~ P1_ID, #cases 3,6
	                           index %in% G1 & P1 %in% G2 & P2  %in% c(G1,G3) ~ P2_ID, #cases 7,8
	                           index %in% G2 & P1 %in% c(G1,G3) & P2  %in% G2 ~ P1_ID, #cases 12,15
	                           index %in% G2 & P1 %in% G2 & P2  %in% c(G1,G3) ~ P2_ID, #cases 16,17
	                           TRUE ~ "Notphased")) %>% 
	mutate(B_col = case_when(B_InhFrom == P1_ID ~ "#E69F00",
	                         B_InhFrom == P2_ID ~ "#39918C",
	                         TRUE ~ "#999999"))


	AD <- as.data.table(AD)
	setnames(AD,colnames(AD),c("index","P1","P2"))

	AD.anno <- AD%>%
	mutate(index_ale_count=stringr::str_count(as.character(index),",|:"),
	       p1_ale_count=stringr::str_count(as.character(P1),",|:"),
	       p2_ale_count=stringr::str_count(as.character(P2),",|:"))%>%
	mutate(index_REF_RD=sapply(index,"[[",1),
	       index_ALT_RD=sapply(index,"[[",2),
	       p1_REF_RD=sapply(P1,"[[",1),
	       p1_ALT_RD=sapply(P1,"[[",2),
	       p2_REF_RD=sapply(P2,"[[",1),
	       p2_ALT_RD=sapply(P2,"[[",2),
	       pr_count=index_ALT_RD+index_REF_RD,
	       pr_ALT_Freq=index_ALT_RD/(index_ALT_RD+index_REF_RD))%>%
	mutate(likelyDN=ifelse(p1_ALT_RD<2&p2_ALT_RD<2&index_ALT_RD>5&p1_REF_RD>10&p2_REF_RD>10&pr_count>10&pr_ALT_Freq>0.2,"TRUE","FALSE"))
	AD.anno <- AD.anno[,c("pr_count","pr_ALT_Freq","likelyDN")]
	AD.anno <- AD.anno %>% 
	mutate(across(where(is.numeric), ~replace(., is.na(.), 0)))

	merged <- cbind(GT.anno ,AD.anno)
	setnames(merged,c("index","P1","P2"),c(PR_ID,P1_ID,P2_ID))
	GenomicRanges::mcols(vcf.gr) <- merged
	return(vcf.gr)
}

if (!is.na(gvcf)){
	err <- FALSE
	if (refver == "hg19"){
		tryCatch({
			gr <- GenomicRanges::GRanges(toString(gsub("chr", "", chr)), IRanges::IRanges(final_from, final_to))
			vcf <- ReadVCF(gvcf, gr, refver)
		}, error = function(e){
			err <- TRUE
		})
		if (err){
			gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(final_from, final_to))
	  	vcf <- ReadVCF(gvcf, gr, refver)
		}
	} else {
	  gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(final_from, final_to))
	  vcf <- ReadVCF(gvcf, gr, refver)
	}
	vcf <- vcf %>% 
		as.data.frame()
	infrom <- vcf$B_InhFrom %>% 
		unique()
	snp_cols <- vector(length = length(infrom))


	for (j in 1:length(snp_cols)){
	snp_cols[j] <- vcf$B_col[which(vcf$B_InhFrom == infrom[j])[1]]
	names(snp_cols)[j] <- infrom[j]
	options(scipen =10000)
	}
	if (!isTRUE(is_trio)){
		snp_cols = c("#999999", "#999999", "#999999")
		p_baf <- ggplot(vcf, aes(x=start,y=pr_ALT_Freq, col = B_InhFrom))+
			geom_point(pch = 20, size = 5)+
			scale_snp+
			style_snp+
			guides(color = guide_legend(override.aes = list(size = 8)))+
			scale_x_continuous(n.breaks=10)+
			ylab("B-allele frequency")+
			scale_color_manual(values = snp_cols)+
			theme(legend.position = "none")+
			coord_cartesian(xlim=c(final_from, final_to), expand=F)
	}	else {
		p_baf <- ggplot(vcf, aes(x=start,y=pr_ALT_Freq, col = B_InhFrom))+
			geom_point(pch = 20, size = 5)+
			scale_snp+
			style_snp+
			guides(color = guide_legend(override.aes = list(size = 8)))+
			scale_x_continuous(n.breaks=10)+
			ylab("B-allele frequency")+
			scale_color_manual(values = snp_cols)+
			coord_cartesian(xlim=c(final_from, final_to), expand=F)
	}
	merged_p <- plot_grid(p0, p1, p, p_baf, p_par2, p_segdup, p_OMIM, nrow=7, rel_heights=c(10,2,10,10,5,3,3))
} else {
	merged_p <- plot_grid(p0, p1, p, p_par2, p_segdup, p_OMIM, nrow=6, rel_heights=c(10,2,10,4,3,3))
}

# Save the plot
ggsave("/tmp/r_ggplot.png", plot = merged_p, dpi=350, height=24, width=45, bg="white")
