#' @title Estimate the Sequencing Depth Size Factors for Peak Statistics Quantification.
#'
#' @description \code{estimateSeqDepth} estimate sequencing depth size factors for each MeRIP-seq samples.
#' Under default setting, the sequencing depth are estimated by the robust estimator defined in package \code{DESeq}.
#' i.e. the median of the ratios to the row geometric means.
#'
#' @details The function takes the input of a \code{\link{SummarizedExomePeak}} object,
#' and it estimates the sequencing depth size factors by the columns of its \link{assay}.
#'
#' @param sep a \code{\link{SummarizedExomePeak}} object.
#' @param from a \code{character} specify the subset of features for sequencing depth estimation, can be one of \code{c("Background", "Modification", "All")}.
#'
#' \describe{
#'  \item{\strong{\code{Background}}}{
#'  The sequencing depths are estimated from the background control regions. This method could make the IP/input LFC estimates become closer to the true modification proportion.
#'  }
#'
#'  \item{\strong{\code{Modification}}}{
#'  The sequencing depths are estimated from the modification peaks/sites.
#'  }
#'
#'  \item{\strong{\code{All}}}{
#'  The sequencing depths are estimated from both the background and the modification regions.
#'  }
#' }
#'
#' Under default setting, the sequencing depth factors are estimated from the background Background regions.
#'
#' @param ... inherited from \code{\link{estimateSizeFactorsForMatrix}}.
#'
#' @examples
#'
#' ### Load the example SummarizedExomPeak object
#' f1 = system.file("extdata", "sep_ex_mod.rds", package="exomePeak2")
#'
#' sep <- readRDS(f1)
#'
#' ### Estimate the sequencing depth size factors
#' sep <- estimateSeqDepth(sep)
#'
#' @seealso \code{\link{normalizeGC}}
#'
#' @return This function will return a \code{\link{SummarizedExomePeak}} object containing newly estimated sequencing depth size factors.
#'
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
#'
#' @aliases estimateSeqDepth
#'
#' @rdname estimateSeqDepth-methods
#'
#' @export
#'
setMethod("estimateSeqDepth",
          "SummarizedExomePeak",
                       function(sep,
                                from = c("Background","Modification","All"),
                                ...){
  from <- match.arg(from)

  if(from == "Background") {
    Background_peaks_indx <- grepl("control", rownames(sep))
    if(sum(Background_peaks_indx) == 0) {
    warning("Cannot find Background peaks, the size factors are estimated using the modification containing peaks.", call. = F,immediate. = T)
    sep$sizeFactor <- estimateSizeFactorsForMatrix(assay( sep ) )
    } else {
    sep$sizeFactor <- estimateSizeFactorsForMatrix(assay( sep[Background_peaks_indx,] ))
    }

  }
  if(from == "Modification"){
    mod_peaks_indx <- grepl("mod", rownames(sep))
    sep$sizeFactor <- estimateSizeFactorsForMatrix(assay( sep[mod_peaks_indx,] ) )
  }

  if(from == "All"){
    sep$sizeFactor <- estimateSizeFactorsForMatrix(assay( sep ) )
  }
  return(sep)
})

