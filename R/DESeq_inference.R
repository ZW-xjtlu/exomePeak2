#' @title Statistical Inference with DESeq on IP over input ratio.
#'
#' @param count_assay a \code{matrix} of read count. Rows of the matrix represent the peaks, and collumns represent the samples.
#' @param design_IP a character vector with values of "IP" and "input", which will indicate the design of MeRIP-Seq experiment.
#' @param p_cutoff a numeric value of the p value cutoff used in DESeq inference.
#' @param p_adj_cutoff a numeric value of the adjusted p value cutoff used in DESeq inference; if provided, values in \code{p_cutoff} will be ignored.
#' @param count_cutoff an integer value indicating the cutoff of the sum of reads count in a window, inference is only performed on the windows with read count bigger than the cutoff. Default value is 10.
#' @param logFC_meth a non negative numeric value of the log2 fold change cutoff used in DESeq inferene for methylated peaks (IP > input).
#' @param logFC_control a non positive numeric value of the log2 fold change cutoff used in DESeq inferene for control peaks (IP < input).
#' @param control_size an integer value indicating the number of randomly selected control peaks returned, by default, all the other peaks after the count cutoff are used as the control.
#'
#' @description \code{DESeq_inference} conduct inference on log2 fold changes of IP over input using the negative binomial model intruduced in DESeq.
#' @return a list of the index for the significant methylated peaks (IP > input) and control peaks (peaks other than methylation peaks).
#'
#'
#'
#' @importFrom DESeq estimateSizeFactors estimateDispersions newCountDataSet nbinomTest
#' @export
DESeq_inference <- function(count_assay,
                            design_IP,
                            p_cutoff = NULL,
                            p_adj_cutoff = 0.05,
                            count_cutoff = 10,
                            logFC_meth = 0,
                            logFC_control = 0,
                            control_size = NULL) {

  stopifnot( !(is.null(p_cutoff) & is.null(p_adj_cutoff)) )

  indx_count <- which( rowSums(count_assay) > count_cutoff )

  cds = newCountDataSet( countData = count_assay[indx_count,],
                         conditions = design_IP )

  cds = estimateSizeFactors( cds ) #Let's just use the default size factors right now.

  cds = estimateDispersions( cds )

  res = nbinomTest( cds, "input", "IP" )

  if( anyNA(res) ){
  res <- as.data.frame(  na.omit(res) )
  }

  if(!is.null(p_adj_cutoff)) {
  stat_sig_indx <- res$padj < p_adj_cutoff
  } else {
  stat_sig_indx <- res$pval < p_cutoff
  }

  sig_peak_meth <- as.numeric( rownames(res)[stat_sig_indx & res$log2FoldChange > logFC_meth] )


  sig_peak_control <- setdiff( as.numeric( rownames(res) ), sig_peak_meth )

  if(!is.null( control_size )) {

  if( length(sig_peak_control) < control_size ) {
    warning("The available number of control peaks is less than defined number.",.Call = FALSE)
  } else {
    sig_peak_control <- sample(sig_peak_control,control_size)
  }
  }

  return(
    list(index_meth = indx_count[as.numeric( sig_peak_meth )],
         index_control = indx_count[as.numeric( sig_peak_control )])
    )
}
