#' @title MeripBamFileList
#'
#' @description An object that summarizes the BAM files used in a MeRIP-Seq experiment.
#'
#' @exportClass MeripBamFileList
#'
#' @export

setClass(
  Class = "MeripBamFileList",
  slots = list(Parameter = "ScanBamParam",
               RandomPrimer = "logical"),
  contains = "BamFileList"
)

#' @title SummarizedExomePeak
#'
#' @description An object that summarizes the methylation sites information,
#' size factors for technical variabilities, and the result of the DESeq2 inference.
#'
#' @import SummarizedExperiment
#' @exportClass SummarizedExomePeak
#'

setClass(
  Class = "SummarizedExomePeak",
  slots = list(DESeq2Results = "data.frame"),
  contains = "RangedSummarizedExperiment"
)

