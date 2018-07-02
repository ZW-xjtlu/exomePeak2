#' @title Quantification and inference of the RNA differential methylation values based on the generalized linear models of negative binomial distribution.
#' @param sep is a summarizedExomePeak object.
#' @param shrinkage_method a character indicating the method for emperical bayes shrinkage, can be one in "apeglm" and "ashr".
#' Please check \code{\link{lfcShrink}} for more information.
#'
#' @param ... inherited arguments from \code{\link{DESeq}}
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
#' @import DESeq2
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
                    shrinkage_method = c("apeglm","ashr"),
                    ...) {

  stopifnot((any(sep$design_Treatment) & any(!sep$design_Treatment)))

  if(is.null(colData( sep )$sizeFactor)) {
    sep <- estimate_size_factors(sep)
  }

  indx_meth <- grepl("meth", rownames( sep ) )

  SE_M <- sep
  SE_M$IPinput = "input"
  SE_M$IPinput[SE_M$design_IP] = "IP"
  SE_M$IPinput = factor(SE_M$IPinput)

  SE_M$Perturbation = "Control"
  SE_M$Perturbation[SE_M$design_Treatment] = "Treatment"
  SE_M$Perturbation = factor(  SE_M$Perturbation )

  if(nrow(GCsizeFactors( sep )) == nrow(sep)) {

    gc_na_indx <- rowSums( is.na(GCsizeFactors(sep)) ) > 0
    #Need to deal with the missing values in GC content size factor.
    Cov = ~ IPinput
    dds = suppressMessages( DESeqDataSet(se = SE_M[(!gc_na_indx) & indx_meth,], design = Cov) )
    normalizationFactors(dds) <- GCsizeFactors(sep)[(!gc_na_indx) & indx_meth,]
    dds$IPinput <- relevel(dds$IPinput, "input")
    dds <- suppressMessages( DESeq(dds) )


    gc_na_indx <- rowSums( is.na(GCsizeFactors(sep)) ) > 0
    Cov = ~ Perturbation + IPinput + Perturbation:IPinput
    dds = suppressMessages( DESeqDataSet(se = SE_M[(!gc_na_indx) & indx_meth,], design = Cov) )
    normalizationFactors(dds) <- GCsizeFactors(sep)[(!gc_na_indx) & indx_meth,]
    dds$IPinput <- relevel(dds$IPinput, "input")
    dds$Perturbation <- relevel(dds$Perturbation, "Control")
    dds <- suppressMessages( DESeq(dds) )


  } else {

    Cov = ~ Perturbation + IPinput + Perturbation:IPinput
    dds = suppressMessages( DESeqDataSet(se = SE_M[indx_meth,], design = Cov) )
    dds$IPinput <- relevel(dds$IPinput, "input")
    dds$Perturbation <- relevel(dds$Perturbation, "Control")
    dds <- suppressMessages( DESeq(dds) )

  }

  #Generation of the DESeq2 report.

  if (nrow(GCsizeFactors( sep )) == nrow(sep)) {

    DS_result  <- lfcShrink(dds=dds, coef=4, type="apeglm")

    quantification_rst <- matrix(NA,nrow = nrow(SE_M[indx_meth,]), ncol = ncol(DS_result))

    colnames(quantification_rst) <- colnames(DS_result)

    quantification_rst <- as.data.frame(quantification_rst)

    quantification_rst[(!gc_na_indx)[indx_meth],] <- as.data.frame( DS_result )

  } else {

   suppressWarnings( quantification_rst <- as.data.frame( lfcShrink(dds=dds, coef=4, type="apeglm") ) )

  }

  rownames(quantification_rst) = rownames(SE_M)[indx_meth]

  DESeq2Results( sep ) = quantification_rst

  return(sep)

})
