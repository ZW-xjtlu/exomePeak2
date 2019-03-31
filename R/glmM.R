#' @title Quantification of RNA modification log2 fold change values based on the generalized linear models of negative binomial distribution.
#'
#' @param sep is a summarizedExomePeak object.
#'
#' @param glm_type a character, that specify the type of generalized linear model used in peak calling; The argument can be one of the "auto", "poisson", "NB", and "DESeq2". Default to be "auto".
#'
#' Under the default setting, the DESeq2 GLM of NB is used on experiments with at least 3 biological replicates for both IP and input samples of the treatment and control.
#' The poisson GLM will be applied in the other cases.
#'
#' @param LFC_shrinkage a character indicating the method for emperical bayes shrinkage on the log2 fold change estimates, can be one in "apeglm", "Gaussian", "ashr", or "none".
#' Please check \code{\link{lfcShrink}} for more information.
#'
#' @param ... Optional arguments passed to \code{\link{DESeq}}
#'
#' @description This function conduct a second round of RNA modification level quantification with generalized linear model estimates of the
#' neative binomial distribution. The RNA modification level is quantified as the log2 IP to input fold change.
#'
#' By default, the modification level estimates will undergoes emperical Bayes shrinkage with a Couchey prior, which is defined in the package \link{apeglm}.
#'
#' @import SummarizedExperiment
#'
#' @import DESeq2
#'
#' @docType methods
#'
#' @name glmM
#'
#' @rdname glmM
#'
#' @export
#'
setMethod("glmM",
          "SummarizedExomePeak",
           function(sep,
                    glm_type = c("auto","poisson", "NB", "DESeq2"),
                    LFC_shrinkage = c("apeglm", "Gaussian", "ashr", "none"),
                    ...) {

  LFC_shrinkage = match.arg(LFC_shrinkage)

  glm_type = match.arg(glm_type)

  if(glm_type == "auto") {
    if( all( table(colData(sep)$design_IP) > 1  ) ) {
      glm_type <- "DESeq2"
    } else {
      glm_type <- "poisson"
    }
  }

  if(glm_type == "poisson") {
    message("Calculating peak statistics with poisson GLM...")
  }

  if(glm_type == "NB") {
    message("Calculating peak statistics with NB GLM...")
  }

  if(glm_type == "DESeq2") {
    message("Calculating peak statistics with DESeq2...")
  }


  stopifnot((any(sep$design_IP) & any(!sep$design_IP)))

  if(is.null(colData( sep )$sizeFactor)) {

    sep <- estimateSeqDepth(sep)

  }

  indx_mod <- grepl("mod", rownames( sep ) )

  SE_M <- sep

  SE_M$IPinput = "input"

  SE_M$IPinput[SE_M$design_IP] = "IP"

  SE_M$IPinput = factor(SE_M$IPinput)

    if(!is.null(GCsizeFactors( sep ))) {

      gc_na_indx <- rowSums( is.na(GCsizeFactors(sep)) ) > 0

      #Need to deal with the missing values in GC content size factor.
      Cov = ~ IPinput

      dds = suppressMessages( DESeqDataSet(se = SE_M[(!gc_na_indx) & indx_mod,], design = Cov) )

      glm_off_sets <- GCsizeFactors(sep)[(!gc_na_indx) & indx_mod,]

      #normalization to make the row geometric means = 0 (since DESeq2 only cares about the difference).
      #and this norm factor is still under the original scale (not log scale glm off set).
      centered_off_sets <- exp(glm_off_sets) / exp(rowMeans(glm_off_sets))

      normalizationFactors(dds) <- centered_off_sets

      rm(glm_off_sets, centered_off_sets)

      dds$IPinput <- relevel(dds$IPinput, "input")


    } else {

      Cov = ~ IPinput

      dds = suppressMessages( DESeqDataSet(se = SE_M[indx_mod,], design = Cov) )

      dds$IPinput <- relevel( dds$IPinput, "input" )

    }

    if(glm_type == "poisson"){
      dispersions(dds) = 0
    }

    if(glm_type == "NB"){
      dds = suppressMessages( estimateDispersions( dds, fitType = "mean" ) )
    }

    if(glm_type == "DESeq2"){
      dds = suppressMessages( estimateDispersions( dds ) )
    }

      dds = suppressMessages( nbinomWaldTest( dds ) )

   #Generation of the DESeq2 report.

    if (!is.null(GCsizeFactors( sep ))) {

      if(LFC_shrinkage == "none") {

        DS_result <- suppressMessages( results( dds, altHypothesis = "greater" ) )

      } else {

        DS_result <- suppressMessages( lfcShrink( dds=dds, res = results( dds, altHypothesis = "greater" ), coef=2, type = LFC_shrinkage  ) )

      }

      quantification_rst <- matrix(NA,nrow = nrow(SE_M[indx_mod,]), ncol = ncol(DS_result))

      colnames(quantification_rst) <- colnames(DS_result)

      quantification_rst <- as.data.frame(quantification_rst)

      quantification_rst[(!gc_na_indx)[indx_mod],] <- as.data.frame( DS_result )

    } else {

      if(LFC_shrinkage == "none") {

      quantification_rst <- suppressMessages( as.data.frame( results( dds, altHypothesis = "greater" ) ) )

      } else {

      quantification_rst <- suppressMessages( as.data.frame( lfcShrink( dds=dds, res = results( dds, altHypothesis = "greater" ), coef=2, type = LFC_shrinkage  ) ) )

      }

    }

  rownames(quantification_rst) = rownames(SE_M)[indx_mod]

  DESeq2Results( sep ) = as.data.frame( quantification_rst )

  return(sep)

})
