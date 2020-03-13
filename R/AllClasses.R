#' @title Maintain and use BAM files in MeRIP-Seq experiment
#'
#' @description An S4 object defined in exomePeak2 that summarizes the BAM files used in a MeRIP-Seq experiment.
#'
#' \code{MeripBamFileList()} provide a convenient format to store and manage
#' the BAM file directories for MeRIP-Seq data set.
#'
#' This class contains \link{BamFileList} from the packge \code{Rsamtools}.
#'
#' @details
#'
#' \strong{Constructors:}
#'
#' \code{MeripBamFileList} can be constructed by \code{scanMeripBAM()}
#'
#' \strong{Accessors:}
#'
#' \code{MeripBamFileList} object share all the accessors with the \link{BamFileList} class,
#' please check it for more information.
#'
#' The frequently used accessors include:
#'
#' \code{metadata()}: Return a list storing the design of MeRIP-Seq experiment.
#'
#' \code{Parameter()}: Access to the BAM FLAG parameters used for BAM file filtering.
#'
#' \code{asMate()}: Return a logical value, TRUE if the BAM file is paired end.
#'
#' It has one additional accessor \code{LibraryType()}
#'
#' \code{LibraryType()} retrieves the strand specificity information of the RNA-Seq library.
#'
#' @examples
#'
#' ### Define BAM File Directories
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
#' f1 = system.file("extdata", "treated_IP1.bam", package="exomePeak2")
#' TREATED_IP_BAM = c(f1)
#' f1 = system.file("extdata", "treated_Input1.bam", package="exomePeak2")
#' TREATED_INPUT_BAM = c(f1)
#'
#' ### For MeRIP-Seq Experiment Without the Treatment Group
#'
#' MeRIP_Seq_Alignment <- scanMeripBAM(
#'  bam_ip = IP_BAM,
#'  bam_input = INPUT_BAM,
#'  paired_end = FALSE
#' )
#'
#' ### For MeRIP-Seq Experiment With the Treatment Group
#'
#' MeRIP_Seq_Alignment <- scanMeripBAM(
#'  bam_ip = IP_BAM,
#'  bam_input = INPUT_BAM,
#'  bam_treated_ip = TREATED_IP_BAM,
#'  bam_treated_input = TREATED_INPUT_BAM,
#'  paired_end = FALSE
#' )
#'
#' LibraryType(MeRIP_Seq_Alignment)
#'
#' Parameter(MeRIP_Seq_Alignment)
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
#'The \code{SummarizedExomePeak} object can be constructed by 3 functions.
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
#' @importFrom methods as is new
#' @importFrom stats pbinom qbinom relevel
#' @importFrom utils capture.output read.table write.csv
#'
#' @export
SummarizedExomePeak <- function(...) new("SummarizedExomePeak", ...)
