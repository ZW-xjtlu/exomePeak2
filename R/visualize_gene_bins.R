#' @title A convenient function to help you quickly visualize the exome bins on a gene.
#' @param gene_id The id of the gene.
#' @param gr_bins A \code{GRanges} object of exomePeak bins, a metadata column indexing for gene id is required.
#' @param grl_exbg A \code{GRangesList} object that define the exon regions of each genes.
#' @return A Gvis diagram.
#'
visualize_gene_bins <- function(gene_id,
                                gr_bins,
                                grl_gene){
require(Gviz)
gr_ex <- grl_gene[[gene_id]]

if(is.null(gr_bins$gene_id)){

gr_b <- subsetByOverlaps(gr_bins, gr_ex)

}else{

gr_b <- gr_bins[gr_bins$gene_id == gene_id]

}

mcols(gr_b) <- NULL

aTrack.stacked <- AnnotationTrack(range = gr_b,
                                  name = paste0( "Bins on Gene ENTREZID: ",
                                                               gene_id),
                                  stacking = "squish",
                                  shape="box")

names(gr_ex) = paste0("exons_", seq_along(gr_ex) )

gr_ex <- keepStandardChromosomes(gr_ex)

grtrack <- AnnotationTrack(gr_ex, genome = "hg19",
                           chromosome = as.vector(seqnames(gr_ex))[1], name = "Gene exons",
                           background.title = "brown",
                           fill="orange",
                           stacking = "dense",
                           shape="box")

plotTracks(list(aTrack.stacked,grtrack), groupAnnotation = "group")

}
