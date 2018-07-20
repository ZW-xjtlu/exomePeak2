#Other experiments for the gcapc approach
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
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

if(!parallel) {
  register(SerialParam())
}


paired <- any( asMates(merip_alignment$BamList) )

exBug <- exons_by_unique_gene(txdb)

org_exons <- exons(txdb)

SE_EX_counts <- summarizeOverlaps(
  features = exBug,
  reads = merip_alignment$BamList,
  param = merip_alignment$Parameter,
  mode = "Union",
  inter.feature = FALSE,
  singleEnd = !paired,
  ignore.strand = !merip_alignment$StrandSpecific,
  fragments = paired
)

saveRDS(SE_EX_counts,"SE_EX_counts.rds")

org_EX_counts <- summarizeOverlaps(
  features = org_exons,
  reads = merip_alignment$BamList,
  param = merip_alignment$Parameter,
  mode = "Union",
  inter.feature = FALSE,
  singleEnd = !paired,
  ignore.strand = !merip_alignment$StrandSpecific,
  fragments = paired
)

saveRDS(org_EX_counts,"SE_exons_counts.rds")


plot_df = reshape2::melt(assay(SE_EX_counts))

plot_df$GC_content = rep(GC_contents$GC_content,10)

p2 <- ggplot(plot_df, aes(x = GC_content,
                          y = value)) +
  geom_point(size = 0.03,
             alpha = 0.1) +
  theme_classic() +
  labs(y = "reads count / region length",
       x = "GC content",
       title = "Reads count scatter plot",
       subtitle = "hg19 exons") +
  scale_y_log10() +
  facet_wrap(~Var2, nrow = 2, ncol = 5) +
  geom_smooth()

#IP / mean(input)

SE_exons_counts <- readRDS("SE_exons_counts.rds")

colData(SE_exons_counts) <- DataFrame( metadata( MeRIP_Seq_Alignment$BamList ) )

GC_contents <- GC_content_over_grl(
  grl = rowRanges(SE_exons_counts),
  bsgenome = Hsapiens,
  fragment_length = 100
)

Plot_df_function <- function(SE,GC,input_cut_off = 10,exon_cut = 400){

mean_input_control <- rowMeans(assay(SE)[,!(SE$design_IP)&!(SE$design_Treatment)]) + 0.01

mean_input_treated <- rowMeans(assay(SE)[,!(SE$design_IP)&(SE$design_Treatment)]) + 0.01

index_keep <- mean_input_control > 10 & mean_input_treated > 10

Ratio_control <- (assay(SE)[index_keep,(SE$design_IP)&!(SE$design_Treatment)])/mean_input_control[index_keep]

Ratio_treated <- (assay(SE)[index_keep,(SE$design_IP)&(SE$design_Treatment)])/mean_input_treated[index_keep]

Ratios <- cbind(Ratio_control,Ratio_treated)

colnames(Ratios) <- c(paste0("control_","rep",1:sum((SE$design_IP)&!(SE$design_Treatment))),
                      paste0("treated_","rep",1:sum((SE$design_IP)&(SE$design_Treatment))))

Plot_df <- reshape2::melt(Ratios)[,-1]

colnames(Plot_df) <- c("IPsamples",
                       "IPinputRatio")

Plot_df$GC_content <- rep(GC$GC_content[index_keep],ncol(Ratios))

Plot_df$Long_exons_index <- rep((GC$Indx_length > exon_cut)[index_keep],ncol(Ratios))

return(Plot_df)
}

plot_df <- Plot_df_function(SE_exons_counts,GC_contents)

p2 <- ggplot(plot_df, aes(x = GC_content,
                          y = IPinputRatio)) +
  geom_point(size = 0.03,
             alpha = 0.1,
             aes(colour = Long_exons_index)) +
  theme_classic() +
  labs(y = "IP reads count / mean input count",
       x = "GC content",
       title = "IP/input ratio scatter plot",
       subtitle = "hg19 exons") +
  scale_y_log10() +
  geom_smooth(aes(group = Long_exons_index,
                  linetype = Long_exons_index),
              colour = "yellow") +
  scale_colour_manual(values = c("blue", "red")) +
  facet_wrap(~IPsamples, ncol = 5)

p3 <- ggplot(plot_df, aes(x = GC_content,
                          y = IPinputRatio)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  theme_classic() +
  labs(y = "IP reads count / mean input count",
       x = "GC content",
       title = "IP/input ratio countor plot",
       subtitle = "hg19 exons") +
  scale_y_log10() +
  facet_wrap(~IPsamples, ncol = 5)

#try reduced exons
GC_contents <- GC_content_over_grl(
  grl = rowRanges(readRDS("SE_EX_counts.rds")),
  bsgenome = Hsapiens,
  fragment_length = 100
)
SE_EX_counts <- readRDS("SE_EX_counts.rds")
colData(SE_EX_counts) <- DataFrame( metadata( MeRIP_Seq_Alignment$BamList ) )

plot_df <- Plot_df_function(SE_EX_counts,GC_contents,exon_cut = 1000)


ggplot(plot_df, aes(x = GC_content,
                    y = IPinputRatio)) +
  geom_point(size = 0.03,
             alpha = 0.1,
             aes(colour = Long_exons_index)) +
  theme_classic() +
  labs(y = "IP reads count / mean input count",
       x = "GC content",
       title = "IP/input ratio scatter plot",
       subtitle = "hg19 exons") +
  scale_y_log10() +
  geom_smooth(aes(group = Long_exons_index,
                  linetype = Long_exons_index),
              colour = "yellow") +
  scale_colour_manual(values = c("blue", "red")) +
  facet_wrap(~IPsamples, ncol = 5)

#Un-informative

SE_exons_counts <- readRDS("SE_exons_counts.rds")

