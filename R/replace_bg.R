#' @title Replace the control regions with user provided background.
#'
#' @param grl A \code{GRangesList} of the merged peaks which the background regions are waiting to be replaced.
#' @param bg A \code{GRanges} or \code{GRangesList} object of the user provided background.
#' @param txdb A \code{TxDb} object that define the transcript annotation.
#' @return A \code{GRangesList} object.
#' The first portion is the exons regions that is not overlapped with \code{annoation}.
#'
#' If the resulting ranges have less number and width compared with what defined in \code{cut_off_width} or \code{cut_off_num},
#' the exon regions of txdb will be returned as the background.
#'
#' The second portion is the restructed user provided annotation with gene id annotated.
#'
#' @import GenomicRanges
#'
#' @importFrom S4Vectors queryHits subjectHits
#'
#'

replace_bg <- function(grl,
                       bg,
                       txdb) {
  grl <- grl[!grepl("control",names(grl))]
  bg_gr <- unlist(bg)
  mcols(bg_gr) <- NULL
  exbyug <- exons_by_unique_gene(txdb)
  bg_gr$gene_id = NA
  fol <- findOverlaps(bg_gr,exbyug)
  bg_gr$gene_id[ queryHits( fol ) ] = names(exbyug)[ subjectHits( fol ) ]
  bg_grl <- split(bg_gr,names(bg_gr))
  names(bg_grl) <- paste0("control_",seq_along(bg_grl))
  return(c(bg_grl,grl))
}
