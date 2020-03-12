#' @title Removing introns from a provided GRanges object.
#' @param gr_bins A \code{GRanges} object of exomePeak bins before the intron removal.
#' @param grl_exbg A \code{GRangesList} object that define the exon regions of each genes.
#' @return A \code{GRangesList} object with the same length of \code{gr_bins}, each list element corresponds to the original GRanges after the removal of introns.
#'
#' @import GenomicRanges
#'
#' @importFrom S4Vectors queryHits subjectHits
#' @keywords internal
#'

remove_introns <- function(gr_bins, grl_exbg){
  #Calculate intronic regions
  Introns_iranges <- gaps(ranges(grl_exbg))
  unlist_ebg <- unlist(grl_exbg)

  seq_lev <- tapply(as.vector( seqnames(unlist_ebg) ), names(unlist_ebg), function(x) x[1] )
  strand_lev <- tapply(as.vector( strand(unlist_ebg) ), names(unlist_ebg), function(x) x[1] )

  #Find the mapping between introns and bins, for only those bins that "contain" introns.
  introns_granges <- GRanges(
    seqnames = rep(seq_lev, elementNROWS(Introns_iranges)),
    ranges = unlist(Introns_iranges),
    strand = rep(strand_lev, elementNROWS(Introns_iranges))
  )

  fol <- findOverlaps(introns_granges,
                      gr_bins,
                      type = "within")

  #Remove all the hits that are inter-genes.
  indx_keep <- names(introns_granges)[queryHits(fol)] == gsub("\\.[0-9]*$","",names(grl_exbg))[gr_bins$transcriptsHits[subjectHits(fol)]]
  fol <- fol[indx_keep,]

  #Split, and re-define the start and ends of those hitted bins.
  indx_Hitted_bins <-  subjectHits(fol)

  bins_contain_introns <- gr_bins[indx_Hitted_bins]
  mcols(bins_contain_introns) <- NULL
  names(bins_contain_introns) <- indx_Hitted_bins

  #For each element within this GRanges, there is going to be an intron / multiple introns.

  introns_each_bins <- introns_granges[queryHits(fol)]
  names(introns_each_bins) <- indx_Hitted_bins

  bins_contain_introns <- c(bins_contain_introns,introns_each_bins)
  bins_contain_introns <- split(bins_contain_introns,names(bins_contain_introns))

  #for some reason, disjoin is slow for huge granges list.

  if(length(bins_contain_introns) == 0) {

    bins_intron_removed <- gr_bins
    mcols(bins_intron_removed) <- NULL
    bins_intron_removed$gene_id <- names(grl_exbg)[gr_bins$transcriptsHits]
    return(bins_intron_removed)

  }else{
    chunk_num = 1e5
    index_start = 1
    for(i in seq_len(ceiling( length(bins_contain_introns)/chunk_num ))) {
      Indx <- index_start: min(i*chunk_num, length(bins_contain_introns))
      bins_contain_introns[Indx] <- disjoin(bins_contain_introns[Indx])
      index_start = i*chunk_num + 1
    }

    #Remove the introns from those GRanges list.
    bins_contain_introns <- unlist(bins_contain_introns)
    bins_contain_introns <- subsetByOverlaps(bins_contain_introns,
                                             introns_granges,
                                             type = "equal",invert = TRUE)
    indx_non_introns <- which( !seq_along(gr_bins) %in% indx_Hitted_bins )

  bins_without_granges <- gr_bins[indx_non_introns]
  mcols(bins_without_granges) <- NULL
  names(bins_without_granges) <- indx_non_introns

  bins_intron_removed <- c(bins_without_granges,bins_contain_introns)

  rm(bins_without_granges)

  rm(bins_contain_introns)

  bins_intron_removed <- bins_intron_removed[order(as.numeric(names(bins_intron_removed)))]

  bins_intron_removed$gene_id <- names(grl_exbg)[gr_bins$transcriptsHits[as.integer( names(bins_intron_removed) )]]

  names(bins_intron_removed) <- names(gr_bins)[as.integer( names(bins_intron_removed) )]

  return(bins_intron_removed)
  }
}
