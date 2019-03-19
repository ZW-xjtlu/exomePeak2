#' @title Estimation of the size factors used in peaks quantification.
#'
#' @description \code{estimateSeqDepth} estimate the sample wised size factors on exomic peaks,
#' by default, it will use the median of the ratios to the row geometric means defined in package DESeq2.
#'
#' @details The function takes the input of a summarizedExomePeak object, and it estimates the size factors based on the count information contained within it.
#'
#' @param sep a summarizedExomePeak object.
#' @param from Determine which group the size factor is estimated from, can be one of the "Control", "Modification", and "Both".
#'
#' By default, the size factors are estimated from the merged control peaks, that correspond to the exonic regions other than the modification containing peaks.
#'
#' @param ... inherited from \code{\link{estimateSizeFactorsForMatrix}}.
#'
#' @return This function will return a summarizedExomePeak object containing newly estimated sample wised size factors.
#'
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
#'
#' @docType methods
#'
#' @name estimateSeqDepth
#'
#' @rdname estimateSeqDepth
#'
#' @export
#'
setMethod("estimateSeqDepth",
          "SummarizedExomePeak",
                       function(sep,
                                from = c("Control","Modification","Both"),
                                ...){
  from <- match.arg(from)

  if(from == "Control") {
    control_peaks_indx <- grepl("control", rownames(sep))
    if(sum(control_peaks_indx) == 0) {
    warning("Cannot find control peaks, the size factors are estimated using the modification containing peaks.", call. = F,immediate. = T)
    sep$sizeFactor <- estimateSizeFactorsForMatrix(assay( sep ) )
    } else {
    sep$sizeFactor <- estimateSizeFactorsForMatrix(assay( sep[control_peaks_indx,] ))
    }

  }
  if(from == "Modification"){
    mod_peaks_indx <- grepl("mod", rownames(sep))
    sep$sizeFactor <- estimateSizeFactorsForMatrix(assay( sep[mod_peaks_indx,] ) )
  }

  if(from == "Both"){
    sep$sizeFactor <- estimateSizeFactorsForMatrix(assay( sep ) )
  }
  return(sep)
})

