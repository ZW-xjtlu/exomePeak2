#' @title Accessor to the slot \code{Parameter} in class \code{MeripBamFileList}.
#'
#' @param x a \code{MeripBamFileList} object.
#'
#' @aliases Parameter
#'
#' @rdname Parameter-methods
#'
#' @examples
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
#' f1 = system.file("extdata", "treated_IP1.bam", package="exomePeak2")
#' TREATED_IP_BAM = c(f1)
#' f1 = system.file("extdata", "treated_Input1.bam", package="exomePeak2")
#' TREATED_INPUT_BAM = c(f1)
#'
#' MeRIP_Seq_Alignment <- scanMeripBAM(
#'   bam_ip = IP_BAM,
#'   bam_input = INPUT_BAM,
#'   paired_end = FALSE
#' )
#'
#' Parameter(MeRIP_Seq_Alignment)
#'
#' @export
#'
#' @return a list for the additional parameters of the MeRIP-seq experiment.
#'
setMethod("Parameter",
          "MeripBamFileList",
          function(x) {
              return(x@Parameter)
          })

#' @title Accessor to the slot \code{LibraryType} in class \code{MeripBamFileList}.
#'
#' @param x a \code{MeripBamFileList} object.
#'
#' @aliases LibraryType
#'
#' @rdname LibraryType-methods
#'
#' @examples
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
#' f1 = system.file("extdata", "treated_IP1.bam", package="exomePeak2")
#' TREATED_IP_BAM = c(f1)
#' f1 = system.file("extdata", "treated_Input1.bam", package="exomePeak2")
#' TREATED_INPUT_BAM = c(f1)
#'
#' MeRIP_Seq_Alignment <- scanMeripBAM(
#'   bam_ip = IP_BAM,
#'   bam_input = INPUT_BAM,
#'   paired_end = FALSE
#' )
#'
#' LibraryType(MeRIP_Seq_Alignment)
#'
#' @return a value for the library type of MeRIP-seq experiment.
#'
#' @export
#'
setMethod("LibraryType",
          "MeripBamFileList",
          function(x) {
              return(x@LibraryType)
          })


#' @title Accessor to the slot \code{GCsizeFactors} in class \code{SummarizedExomePeak}.
#'
#' @param x1 A \code{SummarizedExomePeak} object.
#'
#' @aliases GCsizeFactors
#'
#' @examples
#'
#' f1 = system.file("extdata", "sep_ex_mod.rds", package="exomePeak2")
#'
#' sep <- readRDS(f1)
#'
#' head(GCsizeFactors(sep))
#'
#' @rdname GCsizeFactors-methods
#'
#' @return a data.frame for the GC content size factors of each sample
#'
#' @export
#'
setMethod("GCsizeFactors",
          "SummarizedExomePeak",
          function(x1) {
              return(assays(x1)$GCsizeFactors)
          })

#' @title Accessor to the slot \code{GCsizeFactors} in class \code{SummarizedExomePeak}.
#'
#' @param x2 A \code{SummarizedExomePeak} object.
#' @param value A \code{matrix} object.
#'
#' @aliases GCsizeFactors<-
#'
#' @rdname GCsizeFactors-methods
#'
#' @export
#'
setMethod("GCsizeFactors<-",
          "SummarizedExomePeak",
          function(x2,value) {
              assays(x2,withDimnames=FALSE)$GCsizeFactors <- value
              return(x2)
          })

#' @title Accessor to the slot \code{exomePeak2Results} in class \code{SummarizedExomePeak}.
#'
#' @param x1 A \code{data.frame} object.
#'
#' @aliases exomePeak2Results
#'
#' @rdname exomePeak2Results-methods
#'
#'
#' @examples
#'
#' f1 = system.file("extdata", "sep_ex_mod.rds", package="exomePeak2")
#'
#' sep <- readRDS(f1)
#'
#' head(exomePeak2Results(sep))
#'
#' @export
#'
setMethod("exomePeak2Results",
          "SummarizedExomePeak",
          function(x1) {
              return(x1@exomePeak2Results)
          })

#' @title Accessor to the slot \code{exomePeak2Results} in class \code{SummarizedExomePeak}.
#'
#' @param x2 A \code{SummarizedExomePeak} object.
#' @param value a \code{data.frame} object for the DESeq2 Results.
#'
#' @return A \code{data.frame} object for the DESeq2 Results.
#'
#' @aliases exomePeak2Results<-
#'
#' @rdname exomePeak2Results-methods
#'
#' @export
#'
setMethod("exomePeak2Results<-",
          "SummarizedExomePeak",
          function(x2,value) {
              x2@exomePeak2Results <- value
              return(x2)
          })

