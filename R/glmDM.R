#' @title Quantification and inference of the RNA differential methylation values based on the generalized linear models of negative binomial distribution.
#'
#' @param sep is a summarizedExomePeak object.
#'
#' @param glm_type a character, which can be one of the "auto", "poisson", "NB", and "DESeq2". This argument specify the type of generalized linear model used in peak calling; Default to be "auto".
#'
#' Under the default setting, the DESeq2 GLM of NB is used on experiments with at least 3 biological replicates for both IP and input samples.
#' The poisson GLM will be applied otherwise.
#'
#' @param shrinkage_method a character indicating the method for emperical bayes shrinkage, can be one in "apeglm" and "ashr".
#' Please check \code{\link{lfcShrink}} for more information.
#'
#' @param ... Optional arguments passed to \code{\link{DESeq}}
#'
#' @description This function conducts a second round of RNA differential methylation inference based on an interactive generalized linear model of negative binomial distribution.
#'
#' The differential methylation analysis is performed using the following design:
#'
#' log2(Q) = intercept + I(Treatment) + I(IP) + I(IP):I(Treatment).
#'
#' The statistics returned is calculated based on the coefficient estimate of the interactive term I(IP):I(Treated).
#'
#' The resulting RNA differential methyltion level is quantified in form of the log2 Odds ratio; i.e. log2(IP to input ratio in Treatment / IP to input ratio in Control).
#'
#' By default, the final returned log2 Odds ratio estimate will undergoes emperical Bayes shrinkage using a Couchey prior, which is defined in \link{apeglm}.
#'
#' @import SummarizedExperiment
#'
#' @import DESeq2
#'
#' @docType methods
#'
#' @name glmDM
#'
#' @rdname glmDM
#'
#' @export

setMethod("glmDM",
          "SummarizedExomePeak",
           function(sep,
                    glm_type = c("auto","poisson", "NB", "DESeq2"),
                    shrinkage_method = c("apeglm","ashr","none"),
                    ...) {

  shrinkage_method = match.arg(shrinkage_method)

  glm_type = match.arg(glm_type)

  stopifnot( ( any(sep$design_Treatment) & any(!sep$design_Treatment) ) )

  if(glm_type == "auto") {
    if( all( table(colData(sep)$design_IP,colData(sep)$design_Treatment) > 3 ) ) {
      glm_type <- "DESeq2"
    } else {
      glm_type <- "poisson"
    }
  }

  if(glm_type == "poisson") {
    message("differential modification analysis with poisson GLM")
  }

  if(glm_type == "NB") {
    message("differential modification analysis with NB GLM")
  }

  if(glm_type == "DESeq2") {
    message("differential modification analysis with DESeq2 NB GLM")
  }

  if(is.null(colData( sep )$sizeFactor)) {

    sep <- estimateSeqDepth(sep)

  }

  indx_meth <- grepl("meth", rownames( sep ) )

  SE_M <- sep

  SE_M$IPinput = "input"

  SE_M$IPinput[SE_M$design_IP] = "IP"

  SE_M$IPinput = factor(SE_M$IPinput)

  SE_M$Perturbation = "Control"

  SE_M$Perturbation[SE_M$design_Treatment] = "Treatment"

  SE_M$Perturbation = factor( SE_M$Perturbation )

  if(!is.null(GCsizeFactors( sep ))) {

    gc_na_indx <- rowSums( is.na(GCsizeFactors(sep)) ) > 0

    Cov = ~ Perturbation + IPinput + Perturbation:IPinput

    dds = suppressMessages( DESeqDataSet(se = SE_M[(!gc_na_indx) & indx_meth,], design = Cov) )

    glm_off_sets <- GCsizeFactors(sep)[(!gc_na_indx) & indx_meth,]

    #Normalization to make the row geometric means = 0 (since DESeq2 only cares about the difference)
    #and this norm factor is still under the original scale (not log scale glm off set).

    centered_off_sets <- exp(glm_off_sets) / exp(rowMeans(glm_off_sets))

    normalizationFactors(dds) <- centered_off_sets

    rm(glm_off_sets,centered_off_sets)

  } else {

    Cov = ~ Perturbation + IPinput + Perturbation:IPinput

    dds = suppressMessages( DESeqDataSet(se = SE_M[indx_meth,], design = Cov) )

  }

  dds$IPinput <- relevel(dds$IPinput, "input")

  dds$Perturbation <- relevel(dds$Perturbation, "Control")

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

  #Generation of the DESeq2 report.

  if (!is.null(GCsizeFactors( sep ))) {

   if(shrinkage_method == "none") {

     DS_result <- results(dds)

   } else {

     DS_result  <- lfcShrink(dds = dds, coef = 4, type = shrinkage_method)

   }

    quantification_rst <- matrix( NA, nrow = nrow(SE_M[indx_meth,]), ncol = ncol(DS_result) )

    colnames(quantification_rst) <- colnames(DS_result)

    quantification_rst <- as.data.frame(quantification_rst)

    quantification_rst[(!gc_na_indx)[indx_meth],] <- as.data.frame( DS_result )

  } else {

    if(shrinkage_method == "none") {

    quantification_rst  <- as.data.frame( results(dds) )

    } else {

    suppressWarnings( quantification_rst <- as.data.frame( lfcShrink(dds=dds, coef=4, type = shrinkage_method) ) )

    }

  }

  rownames( quantification_rst ) = rownames( SE_M )[indx_meth]

  DESeq2Results( sep ) = as.data.frame( quantification_rst )

  return( sep )

})
