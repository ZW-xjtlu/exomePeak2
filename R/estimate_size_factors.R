#' @title Estimation of the size factors used in peaks quantification.
#'
#' @description \code{estimate_size_factors} estimate the sample wised size factors on exomic peaks,
#' by default, it will use the median of the ratios to the row geometric means defined in package DESeq2.
#'
#' @details The function takes the input of a summarizedExomePeak object, and it estimates the size factors based on the count information contained within it.
#'
#' @param design_ip an optional logical vector indicating the design of IP and input, with TRUE represents for IP.
#' @param from Determine which group the size factor is estimated from, can be one of the "Control", "Methylation", and "Both".
#'
#' By default, the size factors are estimated from the merged control peaks, that correspond to the exonic regions other than the methylated peaks.
#'
#' @param ... inherited from \code{\link{estimateSizeFactorsForMatrix}}.
#'
#' @return This function will return a summarizedExomePeak object containing newly estimated sample wised size factors.
#'
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
#' @export
estimate_size_factors <- function(sep,
                                  from = c("Control","Methylation","Both"),
                                  ...){
  from <- match.arg(from)

  if(from == "Control") {
    control_peaks_indx <- grepl("control", rownames(sep$SE))
    if(sum(control_peaks_indx) == 0) {
    warning("cannot find control peaks, the size factors are estimated using the methylation peaks.", call. = F,immediate. = T)
    sep$SE$sizeFactor <- estimateSizeFactorsForMatrix(assay( sep$SE ) )
    } else {
    sep$SE$sizeFactor <- estimateSizeFactorsForMatrix(assay( sep$SE[control_peaks_indx,] ))
    }

  }
  if(from == "Methylation"){
    meth_peaks_indx <- grepl("meth", rownames(sep$SE))
    sep$SE$sizeFactor <- estimateSizeFactorsForMatrix(assay( sep$SE[meth_peaks_indx,] ) )
  }

  if(from == "Both"){
    sep$SE$sizeFactor <- estimateSizeFactorsForMatrix(assay( sep$SE ) )
  }
  return(sep)
 }

