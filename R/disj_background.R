#' @title Find the background of the user provided mod_gration.
#'
#' @param mod_gr A \code{GRanges} object of user provided mod_gration (names are neccessary for the index of the splitting).
#' @param txdb A \code{TxDb} object that define the transcript mod_gration.
#' @param cut_off_num  A non-negative integer indicate the leat total number of the disjoint exons used as the background; Default 2000.
#' @param m6Aseq_background A logical value, TRUE if the region of 5'UTR and long exons of the transcripts should be dropped in control region; Default TRUE.
#' @param distance_5p A numeric value of the length of the transcript starting region; default 200.
#' @param rename_meth Whether to rename the returned methylation sites, default = FALSE.
#' @return A \code{GRangesList} object.
#' The first portion is the exons regions that is not overlapped with \code{annoation}.
#'
#' If the resulting ranges have less number and width compared with what defined in \code{cut_off_num},
#' the exon regions of txdb will be returned as the background.
#'
#' The second portion is the restructed user provided mod_gration with gene id mod_grated.
#'
#' @import GenomicRanges
#' @import GenomicFeatures
#'
#' @export

disj_background <- function(mod_gr,
                            txdb,
                            cut_off_num = 2000,
                            background_bins = NULL,
                            background_types = c("mclust", "m6Aseq_prior", "manual", "all"),
                            control_width = 50,
                            rename_meth = FALSE) {

  background_types <- match.arg(background_types)

  exbyug <-
    exons_by_unique_gene(txdb)

  mcols(mod_gr) <- NULL

  mod_gr_tmp <- unlist(mod_gr)

  ######################################################
  #              Prior background of m6A               #
  ######################################################

  if (background_types == "m6Aseq_prior") {

    utr5 <- unlist(fiveUTRsByTranscript(txdb))

    long_exon <- exons(txdb)

    long_exon <- long_exon[width(long_exon) >= 400]

    mcols(utr5) <- NULL

    mcols(long_exon) <- NULL

    mod_gr_tmp <- c(mod_gr_tmp, utr5, long_exon)

    rm(utr5, long_exon)

    disj_ranges <- disjoin(c(unlist(exbyug) , mod_gr_tmp))

    control_ranges <- subsetByOverlaps(disj_ranges,
                                       mod_gr_tmp,
                                       type = "any",
                                       invert = T)

    control_ranges <- reduce(control_ranges)

  } else {

  ######################################################
  #          Background from subset of bins            #
  ######################################################

    bg_bins_tmp <- reduce(unlist(background_bins))

    mcols(bg_bins_tmp) <- NULL

    disj_ranges <- disjoin(c(bg_bins_tmp,
                             mod_gr_tmp))

    rm(bg_bins_tmp)

    control_ranges <- subsetByOverlaps(disj_ranges,
                                       mod_gr_tmp,
                                       type = "any",
                                       invert = T)

    control_ranges <- reduce(control_ranges)
  }


  ######################################################
  #          Check criteria for background             #
  ######################################################

  control_ranges <-
    control_ranges[width(control_ranges) >= control_width]

  if (length(control_ranges) >= cut_off_num) {

    #Annotat the control with gene ids

    control_ranges$gene_id = NA

    fol <- findOverlaps(control_ranges, exbyug)

    control_ranges$gene_id[queryHits(fol)] = names(exbyug)[subjectHits(fol)]

    control_ranges = split(control_ranges,
                           seq_along(control_ranges))

    names(control_ranges) = paste0("control_", names(control_ranges))

  } else {

    #Return background control as total exons

    warning("Not enough exon regions as control, the background used is the total exon.",
            call. = FALSE)

    control_ranges = unlist(exbyug)

    control_ranges$gene_id = names(control_ranges)

    names(control_ranges) <- NULL

    control_ranges = split(control_ranges,
                           seq_along(control_ranges))

    names(control_ranges) = paste0("control_", names(control_ranges))

    seqlengths(control_ranges) = NA

  }

  #organize the mod_grations

  mod_gr$gene_id = NA
  fol <- findOverlaps(mod_gr, exbyug)
  mod_gr$gene_id[queryHits(fol)] = names(exbyug)[subjectHits(fol)]
  split_index <- names(mod_gr)
  names(mod_gr) <- NULL
  mod_grl <- split(mod_gr, split_index)

  if(rename_meth == TRUE)  names(mod_grl) <- seq_along(mod_grl)

  names(mod_grl) <- paste0("meth_", names(mod_grl))

  return(c(mod_grl, control_ranges))

}

