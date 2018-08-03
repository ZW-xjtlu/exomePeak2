#' @title Exome peak calling and quatification of MeRIP-seq datasets.
#'
#' @description This function perform a default exomePeak2 analysis with the following steps:
#'
#'  1. scan MeRIP-seq BAM files: \code{\link{scanMeripBAM}}.
#'
#'  2. call methylation peaks on exons: \code{\link{exomePeakCalling}}.
#'
#'  3. gc correction based on non-modification background regions: \code{\link{GCnormalization}}.
#'
#'  4. inference and quantify modification / differential modification with Negative Binomial GLM:  \code{\link{glmMeth}}; \code{\link{glmDM}}
#'
#'  For complete details on each step, please check the manual pages for the corresponding functions.
#'
#'   call RNA modification peaks and quantify the RNA modification levels on the exon regions.
#' If the treatment samples are provided, differential modification analysis will be conducted to infer the dynamics of RNA modification upon treatment condiftions.
#'
#' @details \code{exomePeak2} call peaks on exonic regions while adjust for GC content biases.
#' A GFF file or a TxDB object should be provided for transcript annotation.
#' Additionally, a BSgenome object is required to perform the GC content adjustment on the methylation peaks.
#'
#' @param merip_bams a \code{MeripBamFileList} object.
#' @param txdb a txdb object, it could be a single character string that can be recognized by \code{\link{makeTxDbFromUCSC}}; Default "hg19".
#' @param fragment_length a positive integer of the expected fragment length in bp; default 100.
#' @param binding_length a positive integer of the antibody binding length in IP samples; default 25.
#' @param step_length a positive integer of the shift size of the sliding window; default is the binding length.
#' @param count_cutoff a non negative integer value of the average reads count per window used in peak calling; default 5.
#' @param p_cutoff a value of the p value cut-off used in peak calling; default NULL.
#' @param p_adj_cutoff a value of the adjusted p value cutoff used in DESeq inference; default 0.05.
#' @param logFC_cutoff a non negative numeric value of the log2 fold change (log2 IP/input) cutoff used in the inferene of peaks.
#' @param width_cutoff a positive integer of the minimum width for the merged peaks; default \code{fragment_length} .
#' @param drop_overlapped_genes a logical indicating whether the bins on overlapping genes are dropped or not; default TRUE.
#' @param parallel a logical indicating whether to use parallel computation, consider this if your computer has more than 16GB RAM.
#' @param mod_annotation a \code{GRanges} object for user provided single based RNA modification annotation. If provided, the peak calling step will be skipped.
#' Reads count will be performed using the provided annotation flanked by length of floor(fragment_length - binding_length/2).
#'
#' The background regions used in this senario will be the disjoint exon regions of the flanked provided sites.
#'
#' @param background a \code{GRanges} or \code{GRangesList} object for user provided background control regions on the genome.
#'
#' @return This function will save the results of modification / differrential modification analysis under the current working directory.
#'
#' @examples
#'
#'
#' @seealso \code{\link{exomePeakCalling}}
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom Rsamtools asMates
#' @import GenomicRanges
#' @import BiocParallel
#' @import SummarizedExperiment
#' @docType methods
#'
#' @name exomePeakCalling
#'
#' @rdname exomePeakCalling
#'
#' @export
#'
exomePeak2 <- function(bam_ip = NULL,
                       bam_input = NULL,
                       bam_treated_ip = NULL,
                       bam_treated_input = NULL,
                       txdb = NULL,
                       bsgenome = NULL,
                       genome_assembly = NA,
                       gene_anno_gff = NULL,
                       paired_end = FALSE,
                       random_primer = TRUE,
                       fragment_length = 100,
                       binding_length = 25,
                       step_length = binding_length,
                       peak_width = fragment_length/2,
                       count_cutoff = 5,
                       p_cutoff = NULL,
                       p_adj_cutoff = 0.05,
                       logFC_cutoff = 0,
                       parallel = FALSE,
                       mod_annotation = NULL,
                       background = NULL,
                       shrinkage_method = c("apeglm","ashr","normal","none"),
                       drop_overlapped_genes = TRUE,
                       gc_correction = TRUE,
                       gc_feature = c("background","all"),
                       export_results = TRUE,
                       export_format = c("tsv","BED","RDS"),
                       table_style = c("bed","granges")
                      ){

shrinkage_method <- match.arg(shrinkage_method)
export_format <- match.arg(export_format)
table_style <- match.arg(table_style)
gc_feature <- match.arg(gc_feature)

if(any(sep$design_Treatment) & shrinkage_method == "normal"){
stop("normal prior is not applicable for differential methylation analysis.")
}

stopifnot(fragment_length > 0)

stopifnot(step_length > 0)

stopifnot(peak_width > 0)

stopifnot(logFC_cutoff >= 0)

stopifnot(count_cutoff >= 0)

if(!is.na(genome_assembly)) {
   if(!is(bsgenome,"BSgenome")) bsgenome = genome_assembly
   if(!is(txdb,"TxDb") & is.null(gene_anno_gff)) txdb = genome_assembly
}

if(!is.null(gene_anno_gff)) {
  txdb <- makeTxDbFromGFF(gene_anno_gff)
} else {
  if (is.null(txdb)) {
    stop("required argument of txdb or gene_anno_gff for transcript annotation.")
  }

  if (!is(txdb, "TxDb")) {
    txdb <- makeTxDbFromUCSC(txdb)
  }
}

if(is.null(bsgenome) & gc_correction) {
  stop("required argument of bsgenome or genome_assemby for GC content correction,
       if you want to skip the GC correction, please let gc_correction = FALSE.")
}

if(gc_correction){
bsgenome <- getBSgenome(bsgenome)
}

#Check the completeness of the genome annotation

merip_bam_lst <- scanMeripBAM(
bam_ip = bam_ip,
bam_input = bam_input,
bam_treated_ip = bam_treated_ip,
bam_treated_input = bam_treated_input,
paired_end = paired_end,
random_primer = random_primer
)

sep <- exomePeakCalling(merip_bams = merip_bam_lst,
                        txdb = txdb,
                        gene_anno_gff = gene_anno_gff,
                        fragment_length = fragment_length,
                        binding_length = binding_length,
                        step_length = binding_length,
                        count_cutoff = count_cutoff,
                        p_cutoff = p_cutoff,
                        p_adj_cutoff = p_adj_cutoff,
                        logFC_cutoff = logFC_cutoff,
                        peak_width = fragment_length/2,
                        drop_overlapped_genes = drop_overlapped_genes,
                        parallel = parallel,
                        mod_annotation = mod_annotation,
                        background = background
)

sep <- estimateSeqDepth(sep)

if(gc_correction) {
  sep <- GCnormalization(sep,
                         bsgenome = bsgenome,
                         txdb = txdb,
                         fragment_length = fragment_length,
                         binding_length = binding_length,
                         drop_overlapped_genes = drop_overlapped_genes,
                         feature = gc_feature)
}

if(any(sep$design_Treatment)){
  sep <- glmDM(sep, shrinkage_method = shrinkage_method)
}else{
  sep <- glmMeth(sep, shrinkage_method = shrinkage_method)
}

if(gc_correction | is(bsgenome, "BSgenome")) {
plotReadsGC(sep = sep,
            save_pdf_prefix = "ep2",
            drop_overlapped_genes = drop_overlapped_genes)

plotBetaGC(sep = sep,
           save_pdf_prefix = "ep2",
           drop_overlapped_genes = drop_overlapped_genes)
}

plotGuitar(sep,
           txdb = txdb,
           save_pdf_prefix = "ep2")

plotExonLength(sep,
               txdb = txdb,
               save_pdf_prefix = "ep2")

if( export_results ) {
exportResults(sep,
              format = export_format,
              inhibit_filter = TRUE,
              table_style = table_style)
}

return(Results(sep,
               table_style = table_style))
}
