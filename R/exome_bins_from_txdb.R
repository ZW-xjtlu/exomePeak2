#' @title extract exome bins for peak calling given a txdb object
#'
#' @param txdb A \code{txdb} object.
#' @param window_size An integer valued number of the width of the sliding windows or bins.
#' @param step_size An integer valued number of the width of the bin steps.
#'
#' @return A \code{GRanges} object of exonic bins with the names corresponding to the indexes of bins.
#' A metadata colomn named gene_id is attached to indicate its gene ID, which is provided by the txdb object.
#' The gene IDs are divided into multiple ones if the gene contains exons that belong to different chromosomes and strands.
#'
#' @import GenomicRanges
#' @import GenomicFeatures
#'
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#'
exome_bins_from_txdb <- function(txdb,
                                 window_size = 25,
                                 step_size = 25) {

  stopifnot(step_size <= window_size)

  #exBygene  <- exonsBy(txdb, by = "gene")

  exBygene  <- exons_by_unique_gene(txdb = txdb)

  tx_widths <- sum(width(exBygene))

  #Try to define the bins start always from the five prime ends of any transcripts / genes.

  bin_nums_on_tx <-
    ceiling(pmax((tx_widths - window_size) / step_size, 1)) + 1 #About 7 million exome bins on hg19.

  strands_tx <- as.vector(strand(unlist(range(exBygene))))

  indx_plus <- strands_tx == "+"

  indx_minus <- strands_tx == "-"

  strands_bins <- rep(strands_tx, bin_nums_on_tx)

  indx_bin_plus <- strands_bins == "+"

  indx_bin_minus <- strands_bins == "-"

  seqnames_bins <- rep(names(tx_widths), bin_nums_on_tx)

  bin_starts_on_tx <- vector("integer", length = sum(bin_nums_on_tx))

  bin_starts_on_tx[indx_bin_plus] <-
    unlist(lapply(bin_nums_on_tx[indx_plus], function(x)
      seq(1, step_size * x, by = step_size)), use.names = FALSE)

  bin_starts_on_tx[indx_bin_minus] <-
    unlist(mapply(
      function(x, y)
        seq(y, y - step_size * (x - 1), by = -1 * step_size),
      bin_nums_on_tx[indx_minus],
      tx_widths[indx_minus]
    ),
    use.names = FALSE) - window_size + 1

  rm(bin_nums_on_tx,
     strands_tx,
     indx_plus,
     indx_minus,
     indx_bin_plus,
     indx_bin_minus)

  bins_on_tx <- GRanges(
    seqnames = seqnames_bins,
    ranges = IRanges(start = bin_starts_on_tx,
                     width = window_size),
    strand = strands_bins
  )

  #Trim over-hanging ends
  tx_widths <- sum(width(exBygene))

  suppressWarnings(seqlengths(bins_on_tx) <-
                     tx_widths[names(seqlengths(bins_on_tx))])

  bins_on_tx <- trim(bins_on_tx)

  bins_on_tx <- bins_on_tx[width(bins_on_tx) >= 10]

  bins_on_genome <-
    suppressWarnings(mapFromTranscripts(bins_on_tx, exBygene))

  names(bins_on_genome) <- seq_along(bins_on_genome)

  rm(bins_on_tx)

  #Removal of introns is time consuming ~ 1min.
  bins_on_genome <-
    remove_introns(bins_on_genome, exBygene) #Ps. some ranges are not preserved, and we need to track them with rownames of the granges.

  # Use the command "visualize_gene_bins("337967",bins_on_genome,exBygene)" to visualize the bins on genome.

  bins_on_genome  <- split(bins_on_genome , names(bins_on_genome))

  bins_on_genome <-
    bins_on_genome[order(as.numeric(names(bins_on_genome)))]

  return(bins_on_genome)

}
