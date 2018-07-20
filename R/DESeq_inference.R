#' @title Statistical Inference with DESeq on IP over input ratio.
#'
#' @param count_assay a \code{matrix} of read count. Rows of the matrix represent the peaks, and the collumns represent the samples.
#' @param design_IP a character vector with values of "IP" and "input", which will indicate the design of MeRIP-Seq experiment.
#' @param p_cutoff a numeric value of the p value cutoff used in DESeq inference.
#' @param p_adj_cutoff a numeric value of the adjusted p value cutoff used in DESeq inference; if provided, values in \code{p_cutoff} will be ignored.
#' @param count_cutoff an integer value indicating the cutoff of the mean of reads count in a row, inference is only performed on the windows with read count bigger than the cutoff. Default value is 10.
#' @param logFC_meth a non negative numeric value of the log2 fold change cutoff used in DESeq inferene for methylated peaks (IP > input).
#' @param min_meth_number a non negative numeric value of the he minimum number of the reported methylated bins.
#' If the bins are filtered less than this number by the p values or effect sizes,
#' more sites will be reported by the order of the p value until it reaches this number; Default to be floor( nrow(count_assay)*0.002 ) .
#'
#' @description \code{DESeq_inference} conduct inference on log2 fold changes of IP over input using the negative binomial model in DESeq.
#' @return a list of the index for the significant methylated peaks (IP > input) and control peaks (peaks other than methylation peaks).
#'
#' @importFrom DESeq estimateSizeFactors estimateDispersions newCountDataSet nbinomTest
DESeq_inference <- function(count_assay,
                            design_IP,
                            p_cutoff = NULL,
                            p_adj_cutoff = 0.05,
                            count_cutoff = 10,
                            logFC_meth = 0,
                            min_meth_number = floor( nrow(count_assay)*0.002 )) {

  stopifnot( !(is.null(p_cutoff) & is.null(p_adj_cutoff)) )

  indx_count <- which( rowMeans(count_assay) > count_cutoff )

  cds = newCountDataSet( countData = count_assay[indx_count,],
                         conditions = design_IP )

  cds = estimateSizeFactors( cds ) #Let's just use the default size factors right now.

  cds = estimateDispersions( cds )

  res = nbinomTest( cds, "input", "IP" )

  if( anyNA(res) ){
  res <- as.data.frame(  na.omit(res) )
  }


  decision_table <- decision_deseq(res = res,
                                   log2FC_cut = logFC_meth,
                                   p_cut =  p_cutoff,
                                   padj_cut = p_adj_cutoff,
                                   min_mod =  min_meth_number)

  if(is.na(decision_table$Cut_Val_expected)) stop("The number of informative bins are too limited to conduct meaningful peak calling.")

  if(!is.null(p_adj_cutoff)) {
  stat_sig_indx <- res$padj < decision_table$Cut_Val_expected
  } else {
  stat_sig_indx <- res$pval < decision_table$Cut_Val_expected
  }

  sig_peak_meth <- as.numeric( rownames(res)[stat_sig_indx & res$log2FoldChange > logFC_meth] )

  return(
  index_meth = indx_count[as.numeric( sig_peak_meth )]
  )
}
