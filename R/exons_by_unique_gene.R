#' @title Extracting exons by the corresponding genes that are on the same chromosomes and strand.
#' @param txdb A TXDB object.
#' @param drop_overlapped_genes A logical value, TRUE if the overlapped genes are deleted from the output.
#' @return A \code{GRangesList} object, each element in it corresponds to GRanges of the reduced exons of an unique gene,
#' the name corresponds to the original gene with .integer indexed if they have exons on different strands and chromosomes.
#'
#' @importFrom GenomicFeatures exonsBy
#' @import GenomicRanges
#' @export
exons_by_unique_gene <- function(txdb,drop_overlapped_genes = TRUE){

rd_exonsBygene <- reduce( exonsBy(txdb, by = "gene") ,min.gapwidth=0L)

genes_ranges <- range( rd_exonsBygene )

unique_genes <- unlist(genes_ranges)

gene_ids_temp <- names(unique_genes)

freq_gene_ids <- table(gene_ids_temp)[unique(gene_ids_temp)]

indx_dup_genes <- rep(freq_gene_ids,freq_gene_ids) > 1

gene_ids_temp[indx_dup_genes] <- paste0(gene_ids_temp[indx_dup_genes],".", sequence(freq_gene_ids[freq_gene_ids > 1]))

names(unique_genes) <- gene_ids_temp

rd_exons <- unlist(rd_exonsBygene)

fol <- findOverlaps(unique_genes,rd_exons)

#remove the "inter gene mapping" in this fol object
indx_keep <- gsub("\\.[0-9]*$","", names(unique_genes)[queryHits(fol)] ) == names(rd_exons)[subjectHits(fol)]

fol <- fol[indx_keep,]

exonsByUniqueGenes <- reduce( split( rd_exons[subjectHits(fol)] , names(unique_genes)[queryHits(fol)] ) ,min.gapwidth=0L)

#Drop the overlapped genes depend on the decisions
if(drop_overlapped_genes){
  exonsByUniqueGenes <-exonsByUniqueGenes[ countOverlaps(exonsByUniqueGenes,exonsByUniqueGenes) == 1]
}

return(exonsByUniqueGenes)
}
