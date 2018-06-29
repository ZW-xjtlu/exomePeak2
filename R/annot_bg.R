#' @title Find the background of the user provided annotation.
#'
#' @param annot A \code{GRanges} object of user provided annotation (names are neccessary for the index of the splitting).
#' @param txdb A \code{TxDb} object that define the transcript annotation.
#' @param cut_off_width A non-negative integer indicate the least total width of the disjoint exons used as the background; Default 1e5.
#' @param cut_off_num  A non-negative integer indicate the leat total number of the disjoint exons used as the background; Default 2000.
#' @param drop_overlapped_genes A logical indicating whether to discard the overlapping genes; Default TRUE.
#' @return A \code{GRangesList} object.
#' The first portion is the exons regions that is not overlapped with \code{annoation}.
#'
#' If the resulting ranges have less number and width compared with what defined in \code{cut_off_width} or \code{cut_off_num},
#' the exon regions of txdb will be returned as the background.
#'
#' The second portion is the restructed user provided annotation with gene id annotated.
#'
#' @import GenomicRanges
#' @export

annot_bg <- function(annot,
                     txdb,
                     cut_off_width = 1e5,
                     cut_off_num = 2000,
                     drop_overlapped_genes = TRUE) {

  #Calculate exon regions
  exbyug <- exons_by_unique_gene(txdb,drop_overlapped_genes = drop_overlapped_genes)

  mcols(annot) <- NULL

  disj_ranges <- disjoin( c( unlist( exbyug ) , unlist(annot) ) )

  control_ranges <- subsetByOverlaps(
                    disj_ranges,
                    annot,
                    type = "any",
                    invert = T )

  if(length(control_ranges) >= cut_off_num & sum(width(control_ranges)) >= cut_off_width) {

   control_ranges$gene_id = NA

  fol <- findOverlaps( control_ranges, exbyug )

   control_ranges$gene_id[ queryHits( fol ) ] = names(exbyug)[ subjectHits( fol ) ]

  control_ranges$gene_id = names(exbyug)[subjectHits( findOverlaps(control_ranges,exbyug) )]

  control_ranges = split(control_ranges,
                         seq_along(control_ranges))

  names(control_ranges) = paste0("control_",names(control_ranges))

  } else {
    #output directly the exons as the background control

    warning("Not enough exon regions as control, the background used is the total exon.",call.=FALSE)

    control_ranges = unlist(exbyug)

    control_ranges$gene_id = names(control_ranges)

    control_ranges = split(control_ranges,
                           seq_along(control_ranges))

    names(control_ranges) = paste0("control_",names(control_ranges))
  }

  annot$gene_id = NA

  fol <- findOverlaps(annot,exbyug)

  annot$gene_id[ queryHits( fol ) ] = names(exbyug)[ subjectHits( fol ) ]

  names(annot) = paste0("meth_",names(annot))

  annot <- split(annot,names(annot))

  return(c(control_ranges,annot))
}
