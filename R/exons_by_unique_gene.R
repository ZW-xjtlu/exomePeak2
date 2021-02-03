#' @title Extracting exons by the corresponding genes that are on the same chromosomes and strand.
#' @param txdb A TXDB object.
#' @return A \code{GRangesList} object, each element in it corresponds to GRanges of the merged exons of an unique gene,
#' the name corresponds to the original gene with .integer indexed if they have exons on different strands and chromosomes.
#' 
#' The genes that are overlap between each other are merged into one genes, and the exons are also merged into one if they can overlap with each other.
#' @importFrom GenomicFeatures exonsBy
#' @import GenomicRanges
#' @keywords internal
exons_by_unique_gene <- function(txdb) {
  #Remove the genes that are not belong to the same strands or chromosomes
  exbg <- exbg[elementNROWS( range(exbg) ) == 1]
  #Creating gene names for duplicated genes
  fol <- findOverlaps(exbg)
  fol <- fol[ queryHits(fol) != subjectHits(fol) ]
  if(nrow(ol_indx_M) == 0){
  rd_exons <- reduce( unlist(exbg), min.gapwidth=0L )
  split_indx[queryHits(fol)] <- names(exbg)[subjectHits(fol)]
  unique_exons_gene <- reduce( split( rd_exons, split_indx ) )
  return(unique_exons_gene)
  #Group reduced exons into relevant genes
  }
  
}
