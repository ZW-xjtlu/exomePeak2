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
  exbg <- exonsBy(txdb, by = "gene")
  exbg <- exbg[elementNROWS(range(exbg)) == 1]
  fol <- findOverlaps(exbg)
  fol <- fol[queryHits(fol) != subjectHits(fol)]
  ol_indx_M <- as.matrix(fol)
  if (nrow(ol_indx_M) == 0) {
    return(exbg)
  }
  else {
    rm(fol)
    new_gene_names_temp <- names(exbg)
    new_gene_names_list <- split(new_gene_names_temp, seq_along(new_gene_names_temp))
    for (i in seq_len(nrow(ol_indx_M))) {
      temp_i <- ol_indx_M[i, 1]
      new_gene_names_list[[temp_i]] <- c(new_gene_names_list[[temp_i]], 
                                         new_gene_names_temp[ol_indx_M[i, 2]])
    }
    rm(ol_indx_M, temp_i, new_gene_names_temp)
    new_gene_names_list <- lapply(new_gene_names_list, sort)
    new_gene_names <- vapply(new_gene_names_list, function(x) paste(x, 
                                                                    collapse = ","), character(1))
    names(exbg) <- new_gene_names
    rm(new_gene_names, new_gene_names_list)
    rd_exons <- reduce(unlist(exbg), min.gapwidth = 0L)
    fol <- findOverlaps(rd_exons, exbg)
    split_indx <- rep(NA, length(rd_exons))
    split_indx[queryHits(fol)] <- names(exbg)[subjectHits(fol)]
    unique_exons_gene <- split(rd_exons, split_indx)
    unique_exons_gene <- reduce(unique_exons_gene)
    return(unique_exons_gene)
  }
}