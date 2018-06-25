#' @title Exome peak calling on MeRIP-seq datasets while considering the biological variabilities.
#'
#' @description \code{merip_peak_calling} call RNA methylation peaks on exome regions with statistical tests that account for the biological variabilities between samples.
#'
#' @details The function conduct exome level peak calling based on the read alignment results and the transcript annotations.
#'
#' @param merip_alignment a \code{MeRIP_seq_alignment} object.
#' @param genome a single character string describing the genome assembly in UCSC, by default it is "hg19". (not implemented yet)
#' @param gff_dir the path of the GFF or GTF file directory for transcript annotation.
#' @param txdb a txdb object, if provided, it will cover the annotation files defined above.
#' @param window_size an integer value of the window or bin width in bp.
#' @param step_size an integer value of the step or shift width in bp.
#' @param count_cutoff a value of the minimum total reads count per window used in peak calling, default to be 10.
#' @param p_cutoff a value of the p value cut-off used in peak calling, default to be 0.05.
#' @param p_adj_cutoff a value of the adjusted p value cutoff used in DESeq inference; if provided, the value of \code{p_cutoff} will be ignored.
#' @param logFC_cutoff a non negative numeric value of the log2 fold change (log2 IP/input) cutoff used in the inferene of peaks.
#' @param drop_overlapped_genes a logical indicating whether the bins on overlapping genes are dropped or not, default to be TRUE.
#' @param parallel a logical indicating whether to use parallel computation, consider this if your computer has more than 16GB RAM.
#' @param provided_annotation a \code{GRanges} or \code{GRangesList} object for user provided RNA modification annotation. If provided, the peak calling will be skipped.
#'
#' Reads counting will be performed using the provided annotation, while all the annotations are resized into the width defined by \code{window_size} fixing at the center.
#'
#' @return This function will return a \code{SummarizedExomePeak} object storing the ranges and reads counts information of the called peaks.
#'
#' @examples
#'
#'
#' @seealso \code{\link{merip_peak_calling}}
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom Rsamtools asMates
#' @import BiocParallel
#' @import SummarizedExperiment
#' @export

merip_peak_calling <- function(merip_alignment = NULL,
                               genome = "hg19",
                               gff_dir = NULL,
                               txdb = NULL,
                               window_size = 100,
                               step_size = 10,
                               count_cutoff = 10,
                               p_cutoff = NULL,
                               p_adj_cutoff = 0.05,
                               logFC_cutoff = 0,
                               drop_overlapped_genes = TRUE,
                               parallel = FALSE,
                               provided_annotation = NULL
                               ){

  if(is.null(provided_annotation)) {

  if(!is.null(gff_dir)){
  txdb <- makeTxDbFromGFF(gff_dir,format = "auto")
  }

  exome_bins_gr <- exome_bins_from_txdb(txdb = txdb,
                                   window_size = window_size,
                                   step_size = step_size,
                                   drop_overlapped_genes = drop_overlapped_genes)

  paired <- any( asMates(merip_alignment$BamList) )

  if(!parallel) {
  register(SerialParam())
  }

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

  rm(exome_bins_gr) #release ~500MB RAM

  colData(SE_Peak_counts) = DataFrame(metadata(merip_alignment$BamList))

  merged_peaks_grl <- call_peaks_with_DESeq(SE_bins = SE_Peak_counts,
                                                count_cutoff = count_cutoff,
                                                p_cutoff = p_cutoff,
                                                p_adj_cutoff = p_adj_cutoff,
                                                logFC_cutoff = logFC_cutoff,
                                                txdb = txdb,
                                                drop_overlapped_genes = drop_overlapped_genes)

  }else{

    merged_peaks_grl = resize( provided_annotation , window_size, fix = "center")

  } #Later need change into the extended POS / tag count

  if(!parallel){
  register(SerialParam())
  }

  SummarizedExomePeaks <- summarizeOverlaps(
    features = merged_peaks_grl,
    reads = merip_alignment$BamList,
    param = merip_alignment$Parameter,
    mode = "Union",
    inter.feature = FALSE,
    singleEnd = !paired,
    ignore.strand = !merip_alignment$StrandSpecific,
    fragments = paired
  )

  colData(SummarizedExomePeaks) = DataFrame(metadata(merip_alignment$BamList))

  SummarizedExomePeaks <- list(SE = SummarizedExomePeaks,
                               FS_sizeFactor = NULL,
                               DESeq2Result = NULL)

  return(SummarizedExomePeaks)
}

