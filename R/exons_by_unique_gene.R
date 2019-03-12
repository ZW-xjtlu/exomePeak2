#' @title Extracting exons by the corresponding genes that are on the same chromosomes and strand.
#' @param txdb A TXDB object.
#' @return A \code{GRangesList} object, each element in it corresponds to GRanges of the reduced exons of an unique gene,
#' the name corresponds to the original gene with .integer indexed if they have exons on different strands and chromosomes.
#'
#' @importFrom GenomicFeatures exonsBy
#' @import GenomicRanges
#'
#' @export
exons_by_unique_gene <- function(txdb) {

exbg <- exonsBy(txdb, by = "gene")

#remove the genes that are not belong to the same strands
exbg <- exbg[elementNROWS( range(exbg) ) == 1]

#Creating gene names for duplicated genes

fol <- findOverlaps(exbg)

fol <- fol[ queryHits(fol) != subjectHits(fol) ]

ol_indx_M <- as.matrix(fol)

rm(fol)

new_gene_names_temp <- names(exbg)

new_gene_names_list <- split(new_gene_names_temp, seq_along(new_gene_names_temp))

#Merge genes that are mutually overlapping
for(i in 1:nrow(ol_indx_M)){
  temp_i <- ol_indx_M[i,1]
  new_gene_names_list[[temp_i]] <- c(new_gene_names_list[[temp_i]],new_gene_names_temp[ol_indx_M[i,2]])
}

rm(ol_indx_M,temp_i,new_gene_names_temp)

new_gene_names_list <- lapply(new_gene_names_list, sort)

new_gene_names <- sapply(new_gene_names_list, function(x) paste(x,collapse = ","))

names(exbg) <- new_gene_names

rm(new_gene_names,new_gene_names_list)

#Group reduced exons into relavent genes

rd_exons <- reduce( unlist(exbg), min.gapwidth=0L )

fol <- findOverlaps( rd_exons, exbg )

split_indx <- rep(NA, length(rd_exons))

split_indx[queryHits(fol)] <- names(exbg)[subjectHits(fol)]

unique_exons_gene <- split( rd_exons, split_indx )

return(unique_exons_gene)
}
