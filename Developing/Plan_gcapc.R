library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(exomePeak2)
library(GenomicAlignments)
library(BiocParallel)

MeRIP_Seq_Alignment <- scan_merip_bams(
  bam_ip = c("./bam/SRR1182619.bam",
             "./bam/SRR1182621.bam",
             "./bam/SRR1182623.bam"),
  bam_input = c("./bam/SRR1182620.bam",
                "./bam/SRR1182622.bam",
                "./bam/SRR1182624.bam"),
  bam_treated_ip = c("./bam/SRR1182603.bam",
                     "./bam/SRR1182605.bam"),
  bam_treated_input = c("./bam/SRR1182604.bam",
                        "./bam/SRR1182606.bam"),
  paired_end = TRUE
)


merip_alignment = MeRIP_Seq_Alignment
genome = "hg19"
gff_dir = NULL
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
window_size = 100
step_size = 10
count_cutoff = 10
p_cutoff = NULL
p_adj_cutoff = 0.05
logFC_cutoff = 0
drop_overlapped_genes = TRUE
parallel = FALSE
provided_annotation = NULL



paired <- any( asMates(merip_alignment$BamList) )

if(!parallel) {
  register(SerialParam())
}

exome_bins_gr <- exome_bins_from_txdb(txdb = txdb,
                                      window_size = window_size,
                                      step_size = step_size,
                                      drop_overlapped_genes = drop_overlapped_genes)

split_with_sorting <- function(ex_bins_gr){
  ex_bins_grl  <- split(ex_bins_gr,names(ex_bins_gr))
  return( ex_bins_grl[ order( as.numeric( names(ex_bins_grl) ) ) ] )
}

SE_Peak_counts <- summarizeOverlaps(
  features = split_with_sorting(exome_bins_gr),
  reads = merip_alignment$BamList,
  param = merip_alignment$Parameter,
  mode = "Union",
  inter.feature = FALSE,
  singleEnd = !paired,
  ignore.strand = !merip_alignment$StrandSpecific,
  fragments = paired
)

saveRDS(SE_Peak_counts, "Bins_100_10.rds")


exome_bins_gr <- exome_bins_from_txdb(txdb = txdb,
                                      window_size = 1000,
                                      step_size = 1000,
                                      drop_overlapped_genes = drop_overlapped_genes)

split_with_sorting <- function(ex_bins_gr){
  ex_bins_grl  <- split(ex_bins_gr,names(ex_bins_gr))
  return( ex_bins_grl[ order( as.numeric( names(ex_bins_grl) ) ) ] )
}

SE_Peak_counts <- summarizeOverlaps(
  features = split_with_sorting(exome_bins_gr),
  reads = merip_alignment$BamList,
  param = merip_alignment$Parameter,
  mode = "Union",
  inter.feature = FALSE,
  singleEnd = !paired,
  ignore.strand = !merip_alignment$StrandSpecific,
  fragments = paired
)

saveRDS(SE_Peak_counts, "Bins_1000_1000.rds")


exome_bins_gr <- exome_bins_from_txdb(txdb = txdb,
                                      window_size = 10000,
                                      step_size = 10000,
                                      drop_overlapped_genes = drop_overlapped_genes)

split_with_sorting <- function(ex_bins_gr){
  ex_bins_grl  <- split(ex_bins_gr,names(ex_bins_gr))
  return( ex_bins_grl[ order( as.numeric( names(ex_bins_grl) ) ) ] )
}

SE_Peak_counts <- summarizeOverlaps(
  features = split_with_sorting(exome_bins_gr),
  reads = merip_alignment$BamList,
  param = merip_alignment$Parameter,
  mode = "Union",
  inter.feature = FALSE,
  singleEnd = !paired,
  ignore.strand = !merip_alignment$StrandSpecific,
  fragments = paired
)


saveRDS(SE_Peak_counts, "Bins_10k_10k.rds")


plot_mm <- function(SE,col_index,bsgenome,save_pdf_prefix){

GC_contents <- GC_content_over_grl(
    grl = rowRanges(SE),
    bsgenome = bsgenome,
    fragment_length = 100
  )

plot_df = data.frame(
  Reads_count = assay(SE)[,col_index],
  GC_content = GC_contents$GC_content
)

p1 <- ggplot(plot_df, aes(x =  GC_content , y = Reads_count )) +
  geom_point(size = .5,
             alpha = .5) +
  theme_classic() +
  labs(x = "GC contents",
       y = "Read counts",
       title = mtitle,
       subtitle = save_pdf_prefix) +
  xlim(c(0.2,0.85)) + scale_y_log10()

if(!is.null( save_pdf_prefix )){

  suppressMessages( ggsave(
    paste0(save_pdf_prefix,
           "_reads_GC.pdf"),
    p1,
    width = 4.5,
    height = 3) )

}

return(p1)
}


SE <- readRDS("Bins_100_10.rds")
#With unbelievable amount of sample sizes

col_index = 1
bsgenome = Hsapiens
save_pdf_prefix = "IP_rep1"

gr <- unlist(rowRanges(SE))
GC_freq <- letterFrequency(Views(bsgenome,gr),letters="CG")
width_each_group <- sum(width(rowRanges(SE)))
sum_freq <- tapply(GC_freq,names(gr),sum)
GC_contents <- sum_freq/width_each_group

GC_contents <- pmin(GC_contents,1)

plot_df = data.frame(
  Reads_count = assay(SE)[,col_index],
  GC_content = GC_contents
)

plot_df <- plot_df[plot_df$Reads_count >= 10,]

p1 <- ggplot(plot_df, aes(x =  GC_content , y = Reads_count )) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  theme_classic() +
  labs(x = "GC contents",
       y = "Read counts",
       title = "Reads count against GC content",
       subtitle = "IP SRR1182619") +
  xlim(c(0.2,0.85)) + scale_y_log10()


index_sum <- rep(1:ceiling( length(sum_freq)/10),each = 10)[1:length(sum_freq)]

GC_contents <- tapply(sum_freq,index_sum,sum)/tapply(width_each_group,index_sum,sum)

GC_contents <- pmin(GC_contents,1)

plot_df = data.frame(
  Reads_count = tapply(assay(SE)[,7],index_sum,sum),
  GC_content = GC_contents
)

plot_df <- plot_df[plot_df$Reads_count >= 1,]

p1 <- ggplot(plot_df, aes(x =  GC_content , y = Reads_count )) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  theme_classic() +
  labs(x = "GC contents",
       y = "Read counts",
       title = "Reads count against GC content",
       subtitle = "IP SRR1182603") +
  xlim(c(0.2,0.85)) + scale_y_log10() + geom_smooth()

SE <- readRDS("Bins_10k_10k.rds")

plot_df = reshape2::melt(assay(SE))

p2 <- ggplot(plot_df, aes(x = value )) +
  geom_density(aes(fill = Var2),
               linetype = 0,
               alpha = 0.4) +
  theme_classic() +
  labs(x = "reads count",
       y = "density",
       title = "Reads count density plot",
       subtitle = "10k bins on exome") +
  scale_x_log10() + facet_wrap(~Var2, nrow = 2, ncol = 5)

GC_contents <- GC_content_over_grl(
  grl = rowRanges(SE),
  bsgenome = bsgenome,
  fragment_length = 100
)

plot_df$GC_content <- rep(GC_contents,10)


hg19_miCLIP <- readRDS("/Users/zhenwei/Documents/GitHub/Project_X/hg19_miCLIP_ml.rds")

library(BSgenome)

GC_contents <- GC_content_over_grl(hg19_miCLIP,Hsapiens)

plot_df <- data.frame(GC_content = GC_contents$GC_content,
                      Group = hg19_miCLIP$Target)

colnames(plot_df)[1] = "GC_content"

boxplot(GC_content~Group,data = plot_df)


plot_df$Group[plot_df$Group == 1] = "m6A RRACH"
plot_df$Group[plot_df$Group == 0] = "random RRACH"

ggplot(plot_df,
       aes(x = GC_content,
           fill = Group)) + geom_density(alpha = 0.5,
                                           linetype = 0) + theme_classic() +
  labs(title = "density of GC content in miCLIP m6A sites",
        x = "GC content")
