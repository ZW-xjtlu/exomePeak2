#' @title A function to reduce ranges between elements within a GRangesList.
#'
#' @description \code{reduce_inter_grl} redivide and reduce the overlapping GRangesList element.
#'
#' @return GRangesList object that is reduced between its inner elements.
#' The metadata collumn is preserved, however the names is not tracked anymore.
#' @param peaks_grl The GRangesList of the peaks region to be reduced.
#' @param txdb The txdb object used during the generation of the peaks.
#'
#' @import GenomicRanges
#' @import GenomicFeatures
#' @keywords internal
#‘
reduce_peaks <- function(peaks_grl,
                         txdb) {
  exBygene  <- exons_by_unique_gene(
    txdb = txdb
  )
  reduced_peaks_on_genome <- mapFromTranscripts( reduce( mapToTranscripts( unlist(peaks_grl) , exBygene) ), exBygene )
  names(reduced_peaks_on_genome) <- reduced_peaks_on_genome$xHits
  reduced_peaks_on_genome <- remove_introns( reduced_peaks_on_genome, exBygene )
  return(reduced_peaks_on_genome)
}
