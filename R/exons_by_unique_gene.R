#' @title Extracting exons by the corresponding genes that are on the same chromosomes and strand.
#' @param txdb A TXDB object.
#' @return A \code{GRangesList} object, each element in it corresponds to GRanges of the merged exons of an unique gene,
#' the name corresponds to the original gene with .integer indexed if they have exons on different strands and chromosomes.
#' 
#' The genes that are overlap between each other are merged into one genes, and the exons are also merged into one if they can overlap with each other.
#' If a gene can belong to 2 stands, the gene will be divided into "2 genes" with their IDs labeled by "_plusStrand" and "_minusStrand" at their ends.
#' 
#' @importFrom GenomicFeatures exonsBy
#' @import GenomicRanges
#' @keywords internal
exons_by_unique_gene <- function(txdb) {
exbg <- exonsBy(txdb, by = "gene")

#add genes that are not belong to the same strands
indx_both_strands <- elementNROWS( range(exbg) ) != 1

if(any(indx_both_strands)){
exbg_both_gr <- unlist(exbg[indx_double])
strands <- as.character(strand(exbg_both_gr))
strands[strands=="+"] <- "_plusStrand"
strands[strands=="-"] <- "_minusStrand"
exbg <- c(exbg[-1*indx_both_strands], split(exbg_both_gr, paste0(names(exbg_both_gr), strands)))
rm(strands,exbg_both_gr)
}

#Creating gene names for duplicated genes
fol <- findOverlaps(exbg)

fol <- fol[ queryHits(fol) != subjectHits(fol) ]

ol_indx_M <- as.matrix(fol)

if(nrow(ol_indx_M) == 0){

  return(exbg)

} else {

  rm(fol)

  new_gene_names_temp <- names(exbg)

  new_gene_names_list <- split(new_gene_names_temp, seq_along(new_gene_names_temp))

  #Merge genes that are mutually overlapping
  for(i in seq_len(nrow(ol_indx_M))){
    temp_i <- ol_indx_M[i,1]
    new_gene_names_list[[temp_i]] <- c(new_gene_names_list[[temp_i]],new_gene_names_temp[ol_indx_M[i,2]])
  }

  rm(ol_indx_M,temp_i,new_gene_names_temp)

  new_gene_names_list <- lapply(new_gene_names_list, sort)

  new_gene_names <- vapply(new_gene_names_list, function(x) paste(x,collapse = ","), character(1))

  names(exbg) <- new_gene_names

  rm(new_gene_names,new_gene_names_list)

  #Group reduced exons into relevant genes

  rd_exons <- reduce( unlist(exbg), min.gapwidth=0L )

  fol <- findOverlaps( rd_exons, exbg )

  split_indx <- rep(NA, length(rd_exons))

  split_indx[queryHits(fol)] <- names(exbg)[subjectHits(fol)]

  unique_exons_gene <- reduce( split( rd_exons, split_indx ) )

  return(unique_exons_gene)
}
}
