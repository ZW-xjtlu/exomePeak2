#' @title Statistical Inference with DESeq package based on the provided reads count for exomic bins.
#'
#' @description \code{call_peaks_with_DESeq} conduct inference on every exome bins using negative binomial model,
#' the significant bins will be the merged into peaks.
#'
#' @details \code{call_peaks_with_DESeq} will performe exome level peak calling using DESeq model,
#' the significant bins in both directions (IP > input and input > IP) are merged separately.
#' The final result will contains both the methylation peaks and control peaks.
#'
#'
#' @param SE_bins a \code{SummarizedExperiment} object. The meta-data collumn should contain the design information of IP & input and treated & control.
#'
#' @param count_cutoff an integer value indicating the cutoff of the mean of reads count in a row, inference is only performed on the windows with read count bigger than the cutoff. Default value is 10.
#'
#' @param p_cutoff a numeric value of the p value cutoff used in DESeq inference.
#'
#' @param p_adj_cutoff a numeric value of the adjusted p value cutoff used in DESeq inference; if provided, the value of \code{p_cutoff} will be ignored.
#'
#' @param logFC_cutoff a non negative numeric value of the log2 fold change (log2 IP/input) cutoff used in the inferene of peaks.
#'
#' @param txdb the txdb object that is necessary for the calculation of the merge of the peaks.
#'
#' @param drop_overlapped_genes A logical indicating whether the overlapping genes were dropped.
#'
#' @return This function will return a list of \code{GRangesList} object storing peaks for both methylation and control results,
#'
#' @import SummarizedExperiment
#'
call_peaks_with_DESeq <- function(SE_bins,
                                  count_cutoff = 10,
                                  p_cutoff = NULL,
                                  p_adj_cutoff = 0.05,
                                  logFC_cutoff = 0,
                                  txdb,
                                  drop_overlapped_genes = TRUE){

  design_IP <- rep("input",ncol(SE_bins))

  design_IP[colData(SE_bins)$design_IP] <- "IP"

  index_meth = DESeq_inference( count_assay = assay(SE_bins),
                               design_IP = design_IP ,
                               p_cutoff = p_cutoff,
                               p_adj_cutoff = p_adj_cutoff,
                               count_cutoff = count_cutoff,
                               logFC_meth = logFC_cutoff )

  gr_meth <- reduce_peaks( peaks_grl = rowRanges(SE_bins)[index_meth],
                           txdb = txdb,
                           drop_overlapped_genes = drop_overlapped_genes )

  return(
    split(gr_meth, names(gr_meth))
 )

}
