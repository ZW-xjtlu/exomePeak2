#' @title Statistical Inference with DESeq2 on IP over input fold change.
#'
#' @param SE_bins a \code{SummarizedExperiment} of read count. It should contains a colData with column named design_IP,
#'  which is a character vector with values of "IP" and "input". The column helps to index the design of MeRIP-Seq experiment.
#' @param glm_type a character, which can be one of  the "Poisson", "NB", and "DESeq2". This argument specify the type of generalized linear model used in peak calling; Default to be "Poisson".
#' The DESeq2 method is only recommended for high power experiments with more than 3 biological replicates for both IP and input.
#' @param p_cutoff a numeric value of the p value cutoff used in DESeq inference.
#' @param p_adj_cutoff a numeric value of the adjusted p value cutoff used in DESeq2 inference; if provided, values in \code{p_cutoff} will be ignored.
#' @param count_cutoff an integer value indicating the cutoff of the mean of reads count in a row, inference is only performed on the windows with read count bigger than the cutoff. Default value is 10.
#' @param logFC_mod a non negative numeric value of the log2 fold change cutoff used in DESeq inferene for modification containing peaks (IP > input).
#' @param min_mod_number a non negative numeric value of the he minimum number of the reported modification containing bins.
#' If the bins are filtered less than this number by the p values or effect sizes,
#' more sites will be reported by the order of the p value until it reaches this number; Default to be floor( nrow(SE_bins)*0.002 ).
#'
#' @param correct_GC_bg a \code{logical} value of whether to estimate the GC content linear effect on background regions; default \code{= FALSE}.
#'
#' If \code{correct_GC_bg = TRUE}, it may result in a more accurate estimation of the technical effect of GC content for the RNA modifications that are highly biologically related to GC content.
#'
#' @param qtnorm a \code{logical} of whether to perform subset quantile normalization after the GC content linear effect correction； default \code{= TRUE}.
#'
#' Subset quantile normalization will be applied within the IP and input samples seperately to account for the inherent differences between the marginal distributions of IP and input samples.
#'
#' @description \code{GLM_inference} conduct inference on log2 fold changes of IP over input using the negative binomial model in DESeq.
#'
#' @return a list of the index for the significant modified peaks (IP > input) and control peaks (peaks other than modification containing peaks).
#'
#' @importFrom DESeq2 DESeqDataSet estimateDispersions estimateSizeFactors nbinomWaldTest results
#'
#'
GLM_inference <- function(SE_bins,
                          glm_type = c("Poisson", "NB", "DESeq2"),
                          p_cutoff = NULL,
                          p_adj_cutoff = 0.01,
                          count_cutoff = 5,
                          logFC_mod = 0,
                          min_mod_number = floor(nrow(SE_bins) * 0.001),
                          correct_GC_bg = FALSE,
                          qtnorm = TRUE) {

  glm_type <- match.arg(glm_type)

  stopifnot(!(is.null(p_cutoff) & is.null(p_adj_cutoff)))

  indx_count <- which(rowMeans(assay(SE_bins)) > count_cutoff)

  dds = DESeqDataSet(se = SE_bins[indx_count, ],
                     design = ~ design_IP)

  ######################################################
  #               Size factor estimation               #
  ######################################################

  dds$sizeFactor = estimateSizeFactorsForMatrix(assay(dds))

  if (!is.null(rowData(SE_bins)$gc_contents)) {
    indx_IP <- dds$design_IP == "IP"

    message("Estimating GC content correction factors for IP samples...")

    if(correct_GC_bg) {
      subindex = which(rowData(dds)$indx_bg & rowData(dds)$indx_gc_est)
    }else{
      subindex = which(rowData(dds)$indx_gc_est)
    }

    cqnObject_IP <- quiet(
      suppressMessages(
      cqn(
        assay(dds)[, indx_IP],
        lengths = rowData(dds)$region_widths,
        lengthMethod = "fixed",
        x = rowData(dds)$gc_contents,
        sizeFactors = dds$sizeFactor[indx_IP],
        subindex = subindex,
        verbose = FALSE,
        sqn = qtnorm
      )
    )
    )

    message("Estimating GC content correction factors for input samples...")

    cqnObject_input <- quiet(
      suppressMessages(
      cqn(
        assay(dds)[, !indx_IP],
        lengths = rowData(dds)$region_widths,
        lengthMethod = "fixed",
        x = rowData(dds)$gc_contents,
        sizeFactors = dds$sizeFactor[!indx_IP],
        subindex = subindex,
        verbose = FALSE,
        sqn = qtnorm
      )
    )
    )

    rm(subindex)

    glm_off_sets <- matrix(NA, nrow = nrow(dds), ncol = ncol(dds))

    rownames(glm_off_sets) = rownames(dds)

    glm_off_sets[, indx_IP] <- cqnObject_IP$glm.offset

    glm_off_sets[, !indx_IP] <- cqnObject_input$glm.offset

    rm(cqnObject_IP, cqnObject_input, indx_IP)

    #normalization to make the row geometric means = 0 (since DESeq2 only cares about the difference).
    #and this norm factor is still under the original scale (not log scale glm off set).
    centered_off_sets <-
      exp(glm_off_sets) / exp(rowMeans(glm_off_sets))

    normalizationFactors(dds) <- centered_off_sets

    rm(glm_off_sets, centered_off_sets)
  }


  ######################################################
  #             Generalized Linear Model               #
  ######################################################

  if (glm_type == "Poisson") {
    dispersions(dds) = 0
  }

  if (glm_type == "NB") {
    dds = estimateDispersions(dds, fitType = "mean")
  }

  if (glm_type == "DESeq2") {
    dds = suppressMessages( estimateDispersions(dds) )
  }

  dds = suppressMessages( nbinomWaldTest(dds) )

  res <- results(dds, altHypothesis = "greater")

  if (anyNA(res)) {
    res <- as.data.frame(na.omit(res))
  }

  if (is.null(p_cutoff)) {
    sig_indx <-
      res$padj < p_adj_cutoff & res$log2FoldChange > logFC_mod

    if (sum(sig_indx) < min_mod_number) {
      warning(
        "The number of positive bins is too small using DESeq2 under current filter,
        the filter is changed into p value < 0.05 and log2FC > 0, please consider using Poisson GLM.",
        call. = FALSE,
        immediate. = TRUE
      )

      sig_indx <- res$pvalue < 0.05 & res$log2FoldChange > 0

      if (sum(sig_indx) < min_mod_number) {
        stop(
          "The bins are not informative to conduct meaningful peak calling, please check the raw data quality."
        )
      }

    }

  } else {
    sig_indx <- res$pvalue < p_cutoff & res$log2FoldChange > logFC_mod

    if (sum(sig_indx) < min_mod_number) {
      warning(
        'The number of positive bins is smaller than the underlimit using DESeq2 method under current filter,
        the filter is changed into p value < 0.05 & log2FC > 0, please consider set glm_type = "Poisson"',
        call. = FALSE,
        immediate. = TRUE
      )

      sig_indx <- res$pvalue < 0.05 & res$log2FoldChange > 0

      if (sum(sig_indx) < min_mod_number) {
        stop(
          "The bins are not informative to conduct meaningful peak calling, please check the raw data quality."
        )
      }

    }

  }

  sig_peak_mod <- as.numeric(rownames(res)[sig_indx])

  return(sig_peak_mod)

}
