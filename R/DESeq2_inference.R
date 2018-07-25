#' @title Statistical Inference with DESeq2 on IP over input fold change.
#'
#' @param SE_bins a \code{SummarizedExperiment} of read count. It should contains a colData with column named design_IP,
#'  which is a character vector with values of "IP" and "input". The column helps to index the design of MeRIP-Seq experiment.
#' @param p_cutoff a numeric value of the p value cutoff used in DESeq inference.
#' @param p_adj_cutoff a numeric value of the adjusted p value cutoff used in DESeq2 inference; if provided, values in \code{p_cutoff} will be ignored.
#' @param count_cutoff an integer value indicating the cutoff of the mean of reads count in a row, inference is only performed on the windows with read count bigger than the cutoff. Default value is 10.
#' @param logFC_meth a non negative numeric value of the log2 fold change cutoff used in DESeq inferene for methylated peaks (IP > input).
#' @param min_meth_number a non negative numeric value of the he minimum number of the reported methylated bins.
#' If the bins are filtered less than this number by the p values or effect sizes,
#' more sites will be reported by the order of the p value until it reaches this number; Default to be floor( nrow(SE_bins)*0.002 ).
#'
#' @description \code{DESeq2_inference} conduct inference on log2 fold changes of IP over input using the negative binomial model in DESeq.
#'
#' @return a list of the index for the significant methylated peaks (IP > input) and control peaks (peaks other than methylation peaks).
#'
#' @importFrom DESeq2 DESeqDataSet estimateDispersions estimateSizeFactors nbinomWaldTest results
DESeq2_inference <- function(SE_bins,
                             p_cutoff = NULL,
                             p_adj_cutoff = 0.05,
                             count_cutoff = 5,
                             logFC_meth = 0,
                             min_meth_number = floor( nrow(SE_bins) * 0.01 )) {

  stopifnot( !(is.null(p_cutoff) & is.null(p_adj_cutoff)) )

  indx_count <- which( rowMeans(assay(SE_bins)) > count_cutoff )

  dds = DESeqDataSet( se = SE_bins[indx_count,],
                      design = ~ design_IP )

  dds = estimateSizeFactors( dds ) #The default size factor estimation method is used in peak calling.

  dds = estimateDispersions( dds )

  dds = nbinomWaldTest( dds )

  res <- results(dds)

  if( anyNA(res) ){
  res <- as.data.frame(  na.omit(res) )
  }

  decision_table <- decision_deseq2(Inf_RES = res,
                                   log2FC_cut = logFC_meth,
                                   P_cut =  p_cutoff,
                                   Padj_cut = p_adj_cutoff,
                                   Min_mod =  min_meth_number,
                                   Exp_dir = "hyper")

  if(is.na(decision_table$Cut_Val_expected)) stop("The number of informative bins are too limited to conduct meaningful peak calling.")

  if(!is.null(p_adj_cutoff)) {

  stat_sig_indx <- res[[decision_table$Cut_By_expected]] < decision_table$Cut_Val_expected

  } else {

  stat_sig_indx <- res[[decision_table$Cut_By_expected]] < decision_table$Cut_Val_expected

  }

  sig_peak_meth <- as.numeric( rownames(res)[stat_sig_indx & res$log2FoldChange > logFC_meth] )

  return(

  sig_peak_meth

  )

}
