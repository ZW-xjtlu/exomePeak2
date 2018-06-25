#' @title Quantification of RNA methylation beta values based on the generalized linear models of negative binomial distribution.
#' @param sep is a summarizedExomePeak object.
#' @param shrinkage_method a character indicating the method for emperical bayes shrinkage, can be one in "apeglm","normal", and "ashr".
#' Please check \code{\link{lfcShrink}} for more information.
#' @param ... inherited arguments from \code{\link{DESeq}}
#' @description This function conduct a second round of RNA methylation level quantification using generalized linear model estimates of the
#' neative binomial distribution. The RNA methyltion level is quantified as the log2 beta-value (log2 IP to input ratio).
#'
#' By default, the methylation level estimates will undergoes emperical Bayes shrinkage using a Couchey prior, which is defined in the package \package{apeglm}.
#'
#' @import SummarizedExperiment
#' @import DESeq2
#' @export
glm_meth <- function(sep,
                     shrinkage_method = c("apeglm","normal","ashr"),
                     ...) {

  stopifnot((any(sep$SE$design_IP) & any(!sep$SE$design_IP)))

  if(is.null(colData( sep$SE )$sizeFactor)){
    sep <- estimate_size_factors(sep)
  }

  indx_meth <- grepl("meth", rownames( sep$SE ) )

  SE_M <- sep$SE
  SE_M$IPinput = "input"
  SE_M$IPinput[SE_M$design_IP] = "IP"
  SE_M$IPinput = factor(SE_M$IPinput)

    if(!is.null(sep$FS_sizeFactor)) {

      gc_na_indx <- rowSums( is.na(sep$FS_sizeFactor) ) > 0
      #Need to deal with the missing values in GC content size factor.
      Cov = ~ IPinput
      dds = suppressMessages( DESeqDataSet(se = SE_M[(!gc_na_indx) & indx_meth,], design = Cov) )
      normalizationFactors(dds) <- sep$FS_sizeFactor[(!gc_na_indx) & indx_meth,]
      dds$IPinput <- relevel(dds$IPinput, "input")
      dds <- suppressMessages( DESeq(dds) )

    } else {

      Cov = ~ IPinput
      dds = suppressMessages( DESeqDataSet(se = SE_M[indx_meth,], design = Cov) )
      dds$IPinput <- relevel(dds$IPinput, "input")
      dds <- suppressMessages( DESeq(dds) )

    }

   #Generation of the DESeq2 report.

    if (!is.null(sep$FS_sizeFactor)) {

      DS_result <- lfcShrink(dds=dds, coef=2, type="apeglm")

      quantification_rst <- matrix(NA,nrow = nrow(SE_M[indx_meth,]), ncol = ncol(DS_result))

      colnames(quantification_rst) <- colnames(DS_result)

      quantification_rst <- as.data.frame(quantification_rst)

      quantification_rst[(!gc_na_indx)[indx_meth],] <- as.data.frame( DS_result )

    } else {

      quantification_rst <- as.data.frame( lfcShrink(dds=dds, coef=2, type="apeglm") )

    }

  rownames(quantification_rst) = rownames(SE_M)[indx_meth]

  sep$DESeq2Result = quantification_rst

  return(sep)
}
