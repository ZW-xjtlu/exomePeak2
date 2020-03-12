#' @title Convert the txdb object into the full transcript and whole genome types
#'
#' @description This function can convert the txdb object into full transcript and whole genome types.
#'
#' @param txdb a \code{TxDb} object containing the regular transcript annotation.
#' @param type the type of \code{TxDb} object of the output, can be one in c("full_tx","whole_genome").
#' @return a TxDb object that will change the exon region into the full transcript and the whole genome regions.
#'
#' @import GenomicFeatures
#'
#' @import GenomicRanges
#'
#' @importFrom IRanges IRanges
#'
#' @importFrom rtracklayer asGFF
#'
#' @name convertTxDb
#' @keywords internal
#'

convertTxDb <- function(txdb,type = c("full_tx", "whole_genome")){

txdb_gff <- asGFF(txdb)
#Remove the previously defined exon and CDS regions
txdb_gff <- txdb_gff[txdb_gff$type != "exon"]
txdb_gff <- txdb_gff[txdb_gff$type != "CDS"]

if(type == "full_tx") {

gr_tx <- transcripts( txdb )

mcols(gr_tx) <- DataFrame(Parent = paste0("TxID:", gr_tx$tx_id),
                          ID = NA,
                          Name = NA,
                          type = "exon")

txdb_gff <- c(txdb_gff,gr_tx)

rm(gr_tx)

} else {

  if(type == "whole_genome"){

    chromosome_lengths <- seqlengths(txdb_gff)

    if(anyNA(chromosome_lengths)){
      chromosome_lengths <- tapply(end(txdb_gff),
                                   as.character(seqnames(txdb_gff)),
                                   max)
    }

    N <- length(chromosome_lengths)

    chrom_names <- names(chromosome_lengths)

    gr_plus <- GRanges(seqnames = chrom_names,
                       ranges = IRanges(start = rep(1, N ),
                                        width = chromosome_lengths),
                      strand = "+")

    gr_minus <- GRanges(seqnames = chrom_names,
                        ranges = IRanges(start = rep(1, N ),
                                        width = chromosome_lengths),
                        strand = "-")

    txdb_gff <- c(gr_plus,gr_minus,
                  gr_plus,gr_minus,
                  gr_plus,gr_minus)

    rm(gr_plus,chromosome_lengths)

    gnames <- c(paste0( chrom_names, "+" ),
                paste0( chrom_names, "-" ))

    mcols(txdb_gff) <- DataFrame(
      type = rep(c("gene","mRNA","exon"), each = 2*N),
      ID = c(paste0("GeneID:",seq_len(2*N)),
              paste0("TxID:",seq_len(2*N)),
              rep(NA,2*N)),
      Name = c( gnames , gnames ,rep(NA,2*N) ),
      Parent = c(rep(NA,2*N),
                  paste0("GeneID:",seq_len(2*N)),
                  paste0("TxID:",seq_len(2*N)))

    )
  }
}

return(makeTxDbFromGRanges(txdb_gff))
}
