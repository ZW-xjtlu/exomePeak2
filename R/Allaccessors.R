#' @title Accessor to the slot \code{Parameter} in class \code{MeripBamFileList}.
#'
#' @param x a \code{MeripBamFileList} object.
#'
#' @aliases Parameter
#'
#' @rdname Parameter-methods
#'
#' @export
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
#' @rdname GCsizeFactors-methods
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
              assays(x2)$GCsizeFactors <- value
              return(x2)
          })

#' @title Accessor to the slot \code{DESeq2Results} in class \code{SummarizedExomePeak}.
#'
#' @param x1 A \code{data.frame} object.
#'
#' @aliases DESeq2Results
#'
#' @rdname DESeq2Results-methods
#'
#' @export
#'
setMethod("DESeq2Results",
          "SummarizedExomePeak",
          function(x1) {
              return(x1@DESeq2Results)
          })

#' @title Accessor to the slot \code{DESeq2Results} in class \code{SummarizedExomePeak}.
#'
#' @param x2 A \code{SummarizedExomePeak} object.
#' @param value A \code{data.frame} object.
#'
#' @aliases DESeq2Results<-
#'
#' @rdname DESeq2Results-methods
#'
#' @export
#'
setMethod("DESeq2Results<-",
          "SummarizedExomePeak",
          function(x2,value) {
              x2@DESeq2Results <- value
              return(x2)
          })

