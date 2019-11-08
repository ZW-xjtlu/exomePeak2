#' @title MeripBamFileList
#'
#' @description An object that summarizes the BAM files used in a MeRIP-Seq experiment.
#'
#' @exportClass MeripBamFileList
#'
#' @aliases MeripBamFileList
#' @name MeripBamFileList-class
#' @rdname MeripBamFileList-class

setClass(
  Class = "MeripBamFileList",
  slots = list(Parameter = "ScanBamParam",
               LibraryType = "character"),
  contains = "BamFileList"
)

#' @title SummarizedExomePeak
#'
#' @description An S4 object defined in exomePeak2 that summarizes the information of modification peaks/sites,
#' reads counts, size factors, GC contents, and the LFC related statistics.
#'
#' This class contains \link{SummarizedExperiment}.
#'
#' @details
#'
#' \strong{Constructors:}
#'
#'The \code{SummarizedExomePeak} object can be contructed by 3 functions.
#'
#'\enumerate{
#'  \item \link{SummarizedExomePeak}
#'  \item \link{exomePeakCalling}
#'  \item \link{exomePeak2}
#'}
#'
#' \strong{Accessors:}
#'
#' The \code{SummarizedExomePeak} object share all the accessors with the \link{SummarizedExperiment} class.
#'
#' It has 2 additional accessors:
#'
#'\enumerate{
#'  \item \link{GCsizeFactors}
#'  \item \link{DESeq2Results}
#'}
#'
#' @examples
#'
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#'
#' # Construct the SummarizedExomePeak object
#'
#' sep <- exomePeak2(bam_ip = c("IP_rep1.bam",
#'                              "IP_rep2.bam",
#'                              "IP_rep3.bam"),
#'                   bam_input = c("input_rep1.bam",
#'                                 "input_rep2.bam",
#'                                 "input_rep3.bam"),
#'                   txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                   bsgenome = Hsapiens)
#'
#' #Access to the slots in the SummarizedExomePeak object
#'
#' ## Access to reads count
#' assay(sep)
#'
#' ## Access to the sequencing depth size factors and experimental design
#' colData(sep)
#'
#' ## Access to the GC content and feature length information
#' elementMetadata(sep)
#'
#' ## Access to the genomic locations of the modification peaks/sites and the background control regions
#' rowRanges(sep)
#'
#' ## Access to the feature specific size factors
#' GCsizeFactors(sep)
#'
#' ## Access to the statistics on (differential) modification LFC
#' DESeq2Results(sep)
#'
#' @import SummarizedExperiment
#' @name SummarizedExomePeak-class
#' @rdname SummarizedExomePeak-class
#' @aliases SummarizedExomePeak
#' @exportClass SummarizedExomePeak
#'

setClass(
  Class = "SummarizedExomePeak",
  slots = list(DESeq2Results = "data.frame"),
  contains = "RangedSummarizedExperiment"
)

#' Wrapper function SummarizedExomePeak.
#' @param ... arguments passed to \code{new()}.
#' @name SummarizedExomePeak
#' @rdname SummarizedExomePeak-class
#' @export
SummarizedExomePeak <- function(...) new("SummarizedExomePeak", ...)
