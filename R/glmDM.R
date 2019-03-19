#' @title Quantification and inference of the RNA differential modification values based on the generalized linear models of negative binomial distribution.
#'
#' @param sep is a summarizedExomePeak object.
#'
#' @param glm_type a character, which can be one of the "auto", "poisson", "NB", and "DESeq2". This argument specify the type of generalized linear model used in peak calling; Default to be "auto".
#'
#' Under the default setting, the DESeq2 GLM of NB is used on experiments with at least 3 biological replicates for both IP and input samples.
#' The poisson GLM will be applied otherwise.
#'
#' @param LFC_shrinkage a character indicating the method for emperical bayes shrinkage on the log2 fold change estimates, can be one in "apeglm" and "ashr".
#' Please check \code{\link{lfcShrink}} for more information.
#'
#' @param ... Optional arguments passed to \code{\link{DESeq}}
#'
#' @description This function conducts a second round of RNA differential modification inference based on an interactive generalized linear model of negative binomial distribution.
#'
#' The differential modification analysis is performed using the following design:
#'
#' log2(Q) = intercept + I(Treatment) + I(IP) + I(IP):I(Treatment).
#'
#' The statistics returned is calculated based on the coefficient estimate of the interactive term I(IP):I(Treated).
#'
#' The resulting RNA differential modification level is quantified in form of the log2 Odds ratio; i.e. log2(IP to input ratio in Treatment / IP to input ratio in Control).
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
                    LFC_shrinkage = c("apeglm","ashr","none"),
                    ...) {

  LFC_shrinkage = match.arg(LFC_shrinkage)

  glm_type = match.arg(glm_type)

  stopifnot( ( any(sep$design_Treatment) & any(!sep$design_Treatment) ) )

  if(glm_type == "auto") {
    if( all( table(colData(sep)$design_IP,colData(sep)$design_Treatment) > 1 ) ) {
      glm_type <- "DESeq2"
    } else {
      glm_type <- "poisson"
    }
  }

  if(glm_type == "poisson") {
    message("Differential modification analysis with poisson GLM...")
  }

  if(glm_type == "NB") {
    message("Differential modification analysis with NB GLM...")
  }

  if(glm_type == "DESeq2") {
    message("Differential modification analysis with DESeq2...")
  }

  if(is.null(colData( sep )$sizeFactor)) {

    sep <- estimateSeqDepth(sep)

  }

  indx_mod <- grepl("mod", rownames( sep ) )

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

    dds = suppressMessages( DESeqDataSet(se = SE_M[(!gc_na_indx) & indx_mod,], design = Cov) )

    glm_off_sets <- GCsizeFactors(sep)[(!gc_na_indx) & indx_mod,]

    #Normalization to make the row geometric means = 0 (since DESeq2 only cares about the difference)
    #and this norm factor is still under the original scale (not log scale glm off set).

    centered_off_sets <- exp(glm_off_sets) / exp(rowMeans(glm_off_sets))

    normalizationFactors(dds) <- centered_off_sets

    rm(glm_off_sets,centered_off_sets)

  } else {

    Cov = ~ Perturbation + IPinput + Perturbation:IPinput

    dds = suppressMessages( DESeqDataSet(se = SE_M[indx_mod,], design = Cov) )

  }

  dds$IPinput <- relevel(dds$IPinput, "input")

  dds$Perturbation <- relevel(dds$Perturbation, "Control")

  if(glm_type == "poisson"){
    dispersions(dds) = 0
  }

  if(glm_type == "NB"){
    dds =  suppressMessages( estimateDispersions( dds, fitType = "mean" ) )
  }

  if(glm_type == "DESeq2"){
    dds = suppressMessages( estimateDispersions( dds ) )
  }


  dds = suppressMessages( nbinomWaldTest( dds ) )

  #Generation of the DESeq2 report.

  if (!is.null(GCsizeFactors( sep ))) {

   if(LFC_shrinkage == "none") {

     DS_result <- suppressMessages( results(dds) )

   } else {

     DS_result  <- suppressMessages( lfcShrink(dds = dds, coef = 4, type = LFC_shrinkage) )

   }

    quantification_rst <- matrix( NA, nrow = nrow(SE_M[indx_mod,]), ncol = ncol(DS_result) )

    colnames(quantification_rst) <- colnames(DS_result)

    quantification_rst <- as.data.frame(quantification_rst)

    quantification_rst[(!gc_na_indx)[indx_mod],] <- as.data.frame( DS_result )

  } else {

    if(LFC_shrinkage == "none") {

    quantification_rst  <- as.data.frame( results(dds) )

    } else {

    suppressWarnings( quantification_rst <- as.data.frame( lfcShrink(dds=dds, coef=4, type = LFC_shrinkage) ) )

    }

  }

  rownames( quantification_rst ) = rownames( SE_M )[indx_mod]

  DESeq2Results( sep ) = as.data.frame( quantification_rst )

  return( sep )

})
