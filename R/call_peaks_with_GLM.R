#' @title Statistical Inference with DESeq package based on the provided reads count for exomic bins.
#'
#' @description \code{call_peaks_with_GLM} conduct inference on every exome bins using negative binomial model,
#' the significant bins will be the merged into peaks.
#'
#' @details \code{call_peaks_with_GLM} will performe exome level peak calling using DESeq2 model,
#'
#' The significant bins will be merged into methylation peaks.
#'
#' The insignificant bins (pass the row means filtering) will also be merged into control peaks.
#'
#' @param SE_bins a \code{SummarizedExperiment} object. The meta-data collumn should contain the design information of IP/input + treated/control.
#'
#' @param glm_type a character, which can be one of  the "poisson", "NB", and "DESeq2". This argument specify the type of generalized linear model used in peak calling; Default to be "poisson".
#' The DESeq2 method is only recommended for high power experiments with more than 3 biological replicates for both IP and input.
#'
#' @param txdb the txdb object that is necessary for the calculation of the merge of the peaks.
#'
#' @param count_cutoff an integer value indicating the cutoff of the mean of reads count in a row, inference is only performed on the windows with the row average read count bigger than the cutoff. Default value is 5.
#'
#' @param p_cutoff a numeric value of the p value cutoff used in DESeq inference.
#'
#' @param p_adj_cutoff a numeric value of the adjusted p value cutoff used in DESeq2 inference; if provided, the value of \code{p_cutoff} will be ignored; Default = 0.05.
#'
#' @param logFC_cutoff a non negative numeric value of the log2 fold change (log2 IP/input) cutoff used in the inferene of peaks.
#'
#' @param drop_overlapped_genes A logical indicating whether the overlapping genes were dropped.
#'
#' @return This function will return a list of \code{GRangesList} object storing peaks for both methylation and control.
#'
#' @import SummarizedExperiment
#'
#'
call_peaks_with_GLM <- function(SE_bins,
                                glm_type = c("poisson","NB","DESeq2"),
                                txdb,
                                count_cutoff = 5,
                                p_cutoff = NULL,
                                p_adj_cutoff = 0.05,
                                logFC_cutoff = 0,
                                drop_overlapped_genes = TRUE) {

  design_IP_temp <- rep( "input", ncol(SE_bins) )

  design_IP_temp[ colData(SE_bins)$design_IP ] <- "IP"

  SE_bins$design_IP <- factor( design_IP_temp )

  SE_bins$design_IP <- relevel( SE_bins$design_IP, "input" )

  rm(design_IP_temp)

  index_meth = GLM_inference( SE_bins = SE_bins,
                              glm_type = glm_type,
                              p_cutoff = p_cutoff,
                              p_adj_cutoff = p_adj_cutoff,
                              count_cutoff = count_cutoff,
                              logFC_meth = logFC_cutoff )

  gr_meth <- reduce_peaks( peaks_grl = rowRanges(SE_bins)[ as.numeric(rownames(SE_bins)) %in% index_meth ],
                           txdb = txdb,
                           drop_overlapped_genes = drop_overlapped_genes )

  return(
    split(gr_meth, names(gr_meth))
 )

}
