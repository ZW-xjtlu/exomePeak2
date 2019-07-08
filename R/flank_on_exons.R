#'@title Flank ranges on exon coordinate.
#'
#'@description This function provides extention of the ends of the GRangesList on transcript coordinate or exons.
#'
#'@param grl a GRangesList object which is the target of the flanking.
#'@param flank_length the length of the flanking regions (on both left and right).
#'@param txdb the TxDb object for the transcript annotation.
#'@param index_flank whether to store the flanking regions seperately; default TRUE.
#'
#'@return A GRanges object with the extented ends on exons.
#' The flanking regions and the binding regions are indexed by the names of the returned GRangesList.
#' The original input GRangesList will be dropped if it is not mapped to the exonic regions.
#'
#'@import GenomicRanges
#'@import GenomicFeatures
#'
#'@export
#'
flank_on_exons <- function(grl,
                           flank_length,
                           txdb,
                           index_flank = TRUE){

  exBygene  <- exons_by_unique_gene(
    txdb = txdb
  )

  bd_on_tx <- mapToTranscripts(unlist(grl), exBygene)

  #remove names of the inner Granges (so don't contain . in the grangeslist name!)

  names(bd_on_tx) <- gsub("\\..*$","",names(bd_on_tx))

  bd_on_tx <- unlist( range( split(bd_on_tx, names(bd_on_tx)) ) )

  if(index_flank) {

    flank_5p <- flank(bd_on_tx, flank_length, start = TRUE)

    names(flank_5p) <- paste0("f5p_", names(bd_on_tx))

    flank_3p <- flank(bd_on_tx, flank_length, start = FALSE)

    names(flank_3p) = paste0("f3p_", names(bd_on_tx))

    names(bd_on_tx) <- paste0("bd_", names(bd_on_tx))

    bins_on_tx <- c(bd_on_tx,flank_5p,flank_3p)

    rm(bd_on_tx,flank_5p,flank_3p)

  } else {

    bins_on_tx <- bd_on_tx + flank_length

    rm(bd_on_tx)

  }

  #Trim over-hanging ends
  tx_widths <- sum( width(exBygene) )

  suppressWarnings( seqlengths(bins_on_tx) <- tx_widths[names(seqlengths(bins_on_tx))] )

  bins_on_tx <- trim(bins_on_tx)

  bins_on_genome <- suppressWarnings( mapFromTranscripts(bins_on_tx,exBygene) )

  rm(bins_on_tx)

  bins_on_genome <- trim( remove_introns(bins_on_genome,exBygene) )

  return(bins_on_genome)

}
