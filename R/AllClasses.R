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
#'  \item \link{exomePeak2Results}
#'}
#'
#' @examples
#'
#' # Generate the SummarizedExomePeak object by peak calling
#'
#' GENE_ANNO_GTF = system.file("extdata", "example.gtf", package="exomePeak2")
#'
#' f1 = system.file("extdata", "IP1.bam", package="exomePeak2")
#' f2 = system.file("extdata", "IP2.bam", package="exomePeak2")
#' f3 = system.file("extdata", "IP3.bam", package="exomePeak2")
#' f4 = system.file("extdata", "IP4.bam", package="exomePeak2")
#' IP_BAM = c(f1,f2,f3,f4)
#' f1 = system.file("extdata", "Input1.bam", package="exomePeak2")
#' f2 = system.file("extdata", "Input2.bam", package="exomePeak2")
#' f3 = system.file("extdata", "Input3.bam", package="exomePeak2")
#' INPUT_BAM = c(f1,f2,f3)
#'
#' sep <- exomePeak2(bam_ip = IP_BAM,
#'                   bam_input = INPUT_BAM,
#'                   gff_dir = GENE_ANNO_GTF,
#'                   genome = "hg19",
#'                   paired_end = FALSE)
#'
#' #Access to the slots in the SummarizedExomePeak object
#'
#' ## Access to reads count
#' assays(sep)
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
#' exomePeak2Results(sep)
#'
#' @import SummarizedExperiment
#' @name SummarizedExomePeak-class
#' @rdname SummarizedExomePeak-class
#' @aliases SummarizedExomePeak
#' @exportClass SummarizedExomePeak
#'
#' @return SummarizedExomePeak object
#'

setClass(
  Class = "SummarizedExomePeak",
  slots = list(exomePeak2Results = "data.frame"),
  contains = "RangedSummarizedExperiment"
)

#' Wrapper function SummarizedExomePeak.
#' @param ... arguments passed to \code{new()}.
#' @name SummarizedExomePeak
#' @rdname SummarizedExomePeak-class
#' @export
SummarizedExomePeak <- function(...) new("SummarizedExomePeak", ...)
