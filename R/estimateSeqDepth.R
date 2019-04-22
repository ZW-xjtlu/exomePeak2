#' @title Estimate the Sequencing Depth Size Factors for Peak Statistics Quantification.
#'
#' @description \code{estimateSeqDepth} estimate sequencing depth size factors for each MeRIP-seq samples.
#' Under default setting, the sequencing depth are estimated by the robust estimator defined in package DESeq2.
#' i.e. the median of the ratios to the row geometric means.
#'
#' @details The function takes the input of a \code{\link{summarizedExomePeak}} object,
#' and it estimates the sequencing depth size factors by the columns of its \link{assay}.
#'
#' @param sep a \code{\link{summarizedExomePeak}} object.
#' @param from a \code{character} specify the subset of features for sequencing depth estimation, can be one of \code{c("Control", "Modification", "Both")}.
#'
#' \describe{
#'  \item{\strong{\code{Control}}}{
#'  The sequencing depths are estimated from the background control regions. This could make the IP/input LFC estimates become a rescaled version of the real modification proportion.
#'  }
#'
#'  \item{\strong{\code{Modification}}}{
#'  The sequencing depths are estimated from the modification peaks/sites.
#'  }
#'
#'  \item{\strong{\code{Both}}}{
#'  The sequencing depths are estimated from both the control and modification features.
#'  }
#' }
#'
#' Under default setting, the sequencing depth factors are estimated from the background control regions.
#'
#' @param ... inherited from \code{\link{estimateSizeFactorsForMatrix}}.
#'
#' @examples
#'
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#'
#' aln <- scanMeripBAM(
#' bam_ip = c("IP_rep1.bam",
#'            "IP_rep2.bam",
#'            "IP_rep3.bam"),
#' bam_input = c("input_rep1.bam",
#'               "input_rep2.bam",
#'               "input_rep3.bam"),
#' paired_end = TRUE
#' )
#'
#'sep <- exomePeakCalling(merip_bams = aln,
#'                        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                        bsgenome = Hsapiens)
#'
#'sep <- estimateSeqDepth(sep)
#'
#' @seealso \code{\link{normalizeGC}}
#'
#' @return This function will return a \code{\link{summarizedExomePeak}} object containing newly estimated sequencing depth size factors.
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

