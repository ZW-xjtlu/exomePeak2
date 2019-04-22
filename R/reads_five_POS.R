#'@title extract reads five prime POS.
#'@param reads a GAlignmentList object.
#'@import GenomicAlignments
#'@import GenomicRanges
#'
#'
reads_five_POS <- function(reads,
                           width = 1,
                           fix = "start",
                           ...) {
  if (is(reads, "GAlignmentsList")) {
    indx_lst <- rep(seq_along(reads), elementNROWS(reads))
    reads <- unlist(reads)
    reads <- as(reads, "GRanges")
    reads_pos <- resize(reads, width = width, fix = fix, ...)
    reads_pos <- split(reads_pos, indx_lst)
    return(reads_pos)

  } else {
    reads <- as(reads, "GRanges")
    reads_pos <- resize(reads, width = width, fix = fix, ...)
    return(reads_pos)
  }
}
