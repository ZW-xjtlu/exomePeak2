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
#' @param window_size a positive integer value of the window or bin width in bp.
#' @param step_size a positive integer value of the step or shift width in bp.
#' @param count_cutoff a non negative integer value of the minimum total reads count per window used in peak calling, default to be 10.
#' @param p_cutoff a value of the p value cut-off used in peak calling, default to be 0.05.
#' @param p_adj_cutoff a value of the adjusted p value cutoff used in DESeq inference; if provided, the value of \code{p_cutoff} will be ignored.
#' @param logFC_cutoff a non negative numeric value of the log2 fold change (log2 IP/input) cutoff used in the inferene of peaks.
#' @param drop_overlapped_genes a logical indicating whether the bins on overlapping genes are dropped or not, default to be TRUE.
#' @param parallel a logical indicating whether to use parallel computation, consider this if your computer has more than 16GB RAM.
#' @param mod_annotation a \code{GRanges} or \code{GRangesList} object for user provided RNA modification annotation. If provided, the peak calling step will be skipped.
#' Reads counting will be performed using the provided annotation, annotations with size less than \code{window_size} are resized to \code{window_size} fixed at the center.
#'
#' The background regions used in this senario will be the exon ranges other than the provided annotation.
#'
#' @param background a \code{GRanges} or \code{GRangesList} object for user provided background control regions on the genome.
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
                               mod_annotation = NULL,
                               background = NULL
                               ){

  stopifnot(window_size > 0)

  stopifnot(step_size > 0)

  stopifnot(logFC_cutoff >= 0)

  stopifnot(count_cutoff >= 0)

  if(is.null(txdb)) {
  if(!is.null(gff_dir)){
  txdb <- makeTxDbFromGFF(gff_dir, format = "auto")
  } else {
    if(!is.null(genome)) {
      txdb <- makeTxDbFromUCSC(genome)
    } else {
      stop("transcript annotation left undefined with any of the arguments: txdb, genome, or gff_dir.")
    }
  }
  }

  paired <- any( asMates(merip_alignment$BamList) )

  if(is.null(mod_annotation)) {

  exome_bins_gr <- exome_bins_from_txdb(txdb = txdb,
                                   window_size = window_size,
                                   step_size = step_size,
                                   drop_overlapped_genes = drop_overlapped_genes)

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

  } else {

    stopifnot( is( mod_annotation, "GRanges") | is( mod_annotation, "GRangesList") )

    if(

      is( mod_annotation, "GRanges")

    ) {

      mod_annotation = split( mod_annotation , 1:length(mod_annotation) )

    }

    mod_annot_gr <- unlist(mod_annotation)

    index_resize <- width( mod_annot_gr  ) < window_size

    mod_annot_gr[index_resize] <- resize( mod_annot_gr[index_resize] , window_size, fix = "center" )

    merged_peaks_grl <- annot_bg( annot = mod_annot_gr,
                                  txdb = txdb,
                                  cut_off_width = 1e5,
                                  cut_off_num = 2000,
                                  drop_overlapped_genes = drop_overlapped_genes )

  }

  if(background) {
    merged_peaks_grl <- replace_bg(grl = merged_peaks_grl,
                                   bg = background,
                                   txdb = txdb,
                                   drop_overlapped_genes = drop_overlapped_genes )
  }

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

  SummarizedExomePeaks <- list(
                               SE = SummarizedExomePeaks,
                               FS_sizeFactor = NULL,
                               DESeq2Result = NULL
                               )

  return(SummarizedExomePeaks)
}

