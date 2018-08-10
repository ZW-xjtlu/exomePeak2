#' @title Statistical Inference with DESeq2 on IP over input fold change.
#'
#' @param SE_bins a \code{SummarizedExperiment} of read count. It should contains a colData with column named design_IP,
#'  which is a character vector with values of "IP" and "input". The column helps to index the design of MeRIP-Seq experiment.
#' @param glm_type a character, which can be one of  the "poisson", "NB", and "DESeq2". This argument specify the type of generalized linear model used in peak calling; Default to be "poisson".
#' The DESeq2 method is only recommended for high power experiments with more than 3 biological replicates for both IP and input.
#' @param p_cutoff a numeric value of the p value cutoff used in DESeq inference.
#' @param p_adj_cutoff a numeric value of the adjusted p value cutoff used in DESeq2 inference; if provided, values in \code{p_cutoff} will be ignored.
#' @param count_cutoff an integer value indicating the cutoff of the mean of reads count in a row, inference is only performed on the windows with read count bigger than the cutoff. Default value is 10.
#' @param logFC_meth a non negative numeric value of the log2 fold change cutoff used in DESeq inferene for methylated peaks (IP > input).
#' @param min_meth_number a non negative numeric value of the he minimum number of the reported methylated bins.
#' If the bins are filtered less than this number by the p values or effect sizes,
#' more sites will be reported by the order of the p value until it reaches this number; Default to be floor( nrow(SE_bins)*0.002 ).
#'
#' @description \code{GLM_inference} conduct inference on log2 fold changes of IP over input using the negative binomial model in DESeq.
#'
#' @return a list of the index for the significant methylated peaks (IP > input) and control peaks (peaks other than methylation peaks).
#'
#' @importFrom DESeq2 DESeqDataSet estimateDispersions estimateSizeFactors nbinomWaldTest results
#' @importFrom MASS glm.nb
#'
GLM_inference <- function(SE_bins,
                          glm_type = c("poisson","NB","DESeq2"),
                          p_cutoff = NULL,
                          p_adj_cutoff = 0.05,
                          count_cutoff = 5,
                          logFC_meth = 0,
                          min_meth_number = floor( nrow(SE_bins) * 0.005 )) {

  glm_type <- match.arg(glm_type)

  stopifnot( !(is.null(p_cutoff) & is.null(p_adj_cutoff)) )

  indx_count <- which( rowMeans(assay(SE_bins)) > count_cutoff )

  dds = DESeqDataSet( se = SE_bins[indx_count,],
                      design = ~ design_IP )

  dds = estimateSizeFactors( dds ) #The default size factor estimation method is used in peak calling.

  if(nrow(rowData(SE_bins)) == nrow(SE_bins)){

    indx_IP <- dds$design_IP == "IP"

    cqnObject_IP <- suppressMessages( cqn(assay(dds)[,indx_IP],
                                          lengths = rowData( dds )$region_widths,
                                          lengthMethod = "smooth",
                                          x = rowData( dds )$gc_contents,
                                          sizeFactors = dds$sizeFactor[indx_IP],
                                          sqn = TRUE,
                                          verbose = FALSE) )

    cqnObject_input <- suppressMessages( cqn(assay(dds)[,!indx_IP],
                                             lengths = rowData( dds )$region_widths,
                                             lengthMethod = "smooth",
                                             x = rowData( dds )$gc_contents,
                                             sizeFactors = dds$sizeFactor[!indx_IP],
                                             sqn = TRUE,
                                             verbose = FALSE) )

    glm_off_sets <- matrix(NA, nrow = nrow(dds), ncol = ncol(dds))

    rownames(glm_off_sets) = rownames(dds)

    glm_off_sets[,indx_IP] <- cqnObject_IP$glm.offset

    glm_off_sets[,!indx_IP] <- cqnObject_input$glm.offset

    rm(cqnObject_IP, cqnObject_input, indx_IP)

    #normalization to make the row geometric means = 0 (since DESeq2 only cares about the difference).
    #and this norm factor is still under the original scale (not log scale glm off set).
    centered_off_sets <- exp(glm_off_sets) / exp(rowMeans(glm_off_sets))

    normalizationFactors(dds) <- centered_off_sets

    rm(glm_off_sets, centered_off_sets)
  }


  if(glm_type == "poisson"){
  dispersions(dds) = 0
  }

  if(glm_type == "NB"){
  dds = estimateDispersions( dds, fitType = "mean" )
  }

  if(glm_type == "DESeq2"){
    dds = estimateDispersions( dds )
  }

  dds = nbinomWaldTest( dds )

  res <- results(dds, altHypothesis = "greater")

  if( anyNA(res) ){
  res <- as.data.frame(  na.omit(res) )
  }

  if(is.null(p_cutoff)) {

  sig_indx <- res$padj < p_adj_cutoff & res$log2FoldChange > logFC_meth

  if(sum(sig_indx) < min_meth_number) {

    warning("The number of positive bins is too small using DESeq2 under current filter, the filter is changed into p value < 0.05 and log2FC > 0, please consider using poisson GLM.",call. = FALSE, immediate. = TRUE)

  sig_indx <- res$pvalue < 0.05 & res$log2FoldChange > 0

  if(sum(sig_indx) < min_meth_number) {stop("The bins are not informative to conduct meaningful peak calling, please check the raw data quality.")}

  }

  } else {

  sig_indx <- res$pvalue < p_cutoff & res$log2FoldChange > logFC_meth

  if(sum(sig_indx) < min_meth_number) {

    warning('The number of positive bins is smaller than the underlimit using DESeq2 method under current filter
         the filter is changed into p value < 0.05 & log2FC > 0, please consider set glm_type = "poisson"',call. = FALSE, immediate. = TRUE)

    sig_indx <- res$pvalue < 0.05 & res$log2FoldChange > 0

    if(sum(sig_indx) < min_meth_number) {stop("The bins are not informative to conduct meaningful peak calling, please check the raw data quality.")}

  }

  }

  sig_peak_meth <- as.numeric( rownames(res)[sig_indx] )


  return(

  sig_peak_meth

  )

}
