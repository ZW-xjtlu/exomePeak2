#'@title extract reads five prime POS and reverse its strand.
#'@param reads a GAlignmentList object.
#'@param width the width of the five prime POS.
#'@param fix the end of the POS, start means 5', end means 3'.
#'@param ... arguments path to function resize.
#'
#'@import GenomicAlignments
#'@import GenomicRanges
#'
#'@return the resized reads GRanges for reversed 5' POS positions.
#'
reads_five_POS_rev <- function(reads,
                           width = 1,
                           fix = "start",
                           ...) {
  if (is(reads, "GAlignmentsList")) {
    indx_lst <- rep(seq_along(reads), elementNROWS(reads))
    reads <- unlist(reads)
    reads <- as(reads, "GRanges")
    reads_pos <- resize(reads, width = width, fix = fix, ...)
    levels( strand(reads_pos) ) = c("-","+","*")
    reads_pos <- split(reads_pos, indx_lst)
    return(reads_pos)

  } else {
    reads <- as(reads, "GRanges")
    reads_pos <- resize(reads, width = width, fix = fix, ...)
    levels( strand(reads_pos) ) = c("-","+","*")
    return(reads_pos)
  }
}
