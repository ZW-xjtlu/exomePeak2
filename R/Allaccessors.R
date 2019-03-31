#' @title Accessor to the slot \code{Parameter} for the object \code{MeripBamFileList}.
#'
#' @param x a \code{MeripBamFileList} object.
#'
#' @return a \code{ScanBamParam} object.
#'
#' @docType methods
#'
#' @name Parameter
#'
#' @rdname Parameter
#'
#' @export
#'
setMethod("Parameter",
          "MeripBamFileList",
          function(x) {
            return(x@Parameter)
          })

#' @title Accessor to the slot \code{LibraryType} for the object \code{MeripBamFileList}.
#'
#' @param x a \code{MeripBamFileList} object.
#'
#' @return a logical vector.
#'
#' @docType methods
#'
#' @name LibraryType
#'
#' @rdname LibraryType
#'
#' @export
#'
setMethod("LibraryType",
          "MeripBamFileList",
          function(x) {
            return(x@LibraryType)
          })


#' @title Accessor to the slot \code{GCsizeFactors} for the object \code{SummarizedExomePeak}.
#'
#' @param x A \code{SummarizedExomePeak} object.
#'
#' @return A matrix.
#'
#' @docType methods
#'
#' @name GCsizeFactors
#'
#' @rdname GCsizeFactors
#'
#' @export
#'
setMethod("GCsizeFactors",
          "SummarizedExomePeak",
          function(x) {
            return(assays(x)$GCsizeFactors)
          })

#' @title Accessor to the slot \code{GCsizeFactors} for the object \code{SummarizedExomePeak}.
#'
#' @param x A \code{SummarizedExomePeak} object.
#'
#' @return A matrix.
#'
#' @docType methods
#'
#' @name GCsizeFactors
#'
#' @rdname GCsizeFactors
#'
#' @export
#'
setMethod("GCsizeFactors<-",
          "SummarizedExomePeak",
          function(x,value) {
            assays(x)$GCsizeFactors <- value
            return(x)
          })

#' @title Accessor to the slot \code{DESeq2Results} for the object \code{SummarizedExomePeak}.
#'
#' @param x A \code{SummarizedExomePeak} object.
#'
#' @return A data.frame.
#'
#' @docType methods
#'
#' @name DESeq2Results
#'
#' @rdname DESeq2Results
#'
#' @export
#'
setMethod("DESeq2Results",
          "SummarizedExomePeak",
          function(x) {
            return(x@DESeq2Results)
          })

#' @title Accessor to the slot \code{DESeq2Results} for the object \code{SummarizedExomePeak}.
#'
#' @param x A \code{SummarizedExomePeak} object.
#'
#' @return A data.frame.
#'
#' @docType methods
#'
#' @name DESeq2Results
#'
#' @rdname DESeq2Results
#'
#' @export
#'
setMethod("DESeq2Results<-",
          "SummarizedExomePeak",
          function(x,value) {
            x@DESeq2Results <- value
            return(x)
          })

