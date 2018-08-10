#' @title Find the background of the user provided annotation.
#'
#' @param annot A \code{GRanges} object of user provided annotation (names are neccessary for the index of the splitting).
#' @param txdb A \code{TxDb} object that define the transcript annotation.
#' @param cut_off_width A non-negative integer indicate the least total width of the disjoint exons used as the background; Default 1e5.
#' @param cut_off_num  A non-negative integer indicate the leat total number of the disjoint exons used as the background; Default 2000.
#' @param drop_overlapped_genes A logical indicating whether to discard the overlapping genes; Default TRUE.
#' @param drop_5p A logical value, TRUE is the region of the five prime start of the transcripts should be dropped in control region; Default TRUE.
#' @param distance_5p A numeric value of the length of the transcript starting region; default 200.
#' @return A \code{GRangesList} object.
#' The first portion is the exons regions that is not overlapped with \code{annoation}.
#'
#' If the resulting ranges have less number and width compared with what defined in \code{cut_off_width} or \code{cut_off_num},
#' the exon regions of txdb will be returned as the background.
#'
#' The second portion is the restructed user provided annotation with gene id annotated.
#'
#' @import GenomicRanges
#' @import GenomicFeatures

annot_bg <- function(annot,
                     txdb,
                     cut_off_width = 1e5,
                     cut_off_num = 2000,
                     drop_overlapped_genes = TRUE,
                     drop_5p = FALSE,
                     distance_5p = 200,
                     control_width = 50) {

  #Calculate exon regions
  exbyug <- exons_by_unique_gene(txdb, drop_overlapped_genes = drop_overlapped_genes)

  mcols(annot) <- NULL

  if(drop_5p){

    tss_5p <- resize( transcripts(txdb), 1, fix = "start")

    TSS_on_tx <- mapToTranscripts(tss_5p, exbyug)

    TSS_on_tx <- resize(TSS_on_tx, distance_5p, fix = "start")

    tss_5p_extended <- mapFromTranscripts(TSS_on_tx, exbyug)

    mcols(tss_5p_extended) <- NULL

    annot_tmp <- c(unlist(annot), tss_5p_extended)

    rm(tss_5p, TSS_on_tx, tss_5p_extended)

  } else {

    annot_tmp <- unlist(annot)

  }

  disj_ranges <- disjoin( c( unlist( exbyug ) , annot_tmp ) )

  control_ranges <- subsetByOverlaps(
                    disj_ranges,
                    annot_tmp,
                    type = "any",
                    invert = T )

  control_ranges <- control_ranges[width(control_ranges) >= control_width]

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

  #organize the annotations

  annot$gene_id = NA

  fol <- findOverlaps(annot,exbyug)

  annot$gene_id[ queryHits( fol ) ] = names(exbyug)[ subjectHits( fol ) ]

  names(annot) = paste0("meth_",names(annot))

  annot <- split(annot,names(annot))

  return(c(annot,control_ranges))
}
