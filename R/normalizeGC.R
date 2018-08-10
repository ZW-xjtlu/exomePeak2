#'@title Estimation of the GC content normalization factors.
#'@param sep a \code{summarizedExomePeak} object.
#'
#'@param bsgenome a \code{\link{BSgenome}} object for the genome sequence, or it could be the name of the reference genome recognized by \code{\link{getBSgenom}}.
#'
#'@param txdb a \code{\link{TxDb}} object for the transcript annotation, or it could be the name of the reference genome recognized by \code{\link{makeTxDbFromUCSC}}.
#'
#'@param fragment_length the expected fragment length of the sequencing library; Default 100.
#'
#'@param binding_length the expected antibody binding length of IP; Default 25.
#'
#'@param feature the features used in the GC effect estimation, can be "background" and "all".
#'If "all" is choosed, the GC effect function will be estimated using both the methylated and the background regions,
#'this choice will force the resulting modification signals independent of GC content, which could be more robust when the background is incorrectly estimated.
#'However,it is possibably biased if the modification level is biologically dependent on GC content;
#'Default "background".
#'
#'@param qtnorm A logical indicating wheather perform quantile normalization after the GC effect estimation； default TRUE.
#'Quantile normalization will be conducted on IP and input samples seperately to account for the biological difference between the marginal distributions of IP and input.
#'
#'@param effective_gc whether to calculate the weighted GC content by the probability of reads alignment; default FALSE.
#'
#'@param drop_overlapped_genes whether to mask the overlapped genes in gene annotation, this is meaningful because the GC content estimation is conducted only on exons; Default TRUE.
#'
#'@description This function estimates the feature specific size factors in order to correct the GC content artifacts,
#'the GC content biases can contributed to the following sources of variabilities in MeRIP-seq data:
#'
#'1. Technical variablities between replicates.
#'
#'2. Batch effects between different laboratories.
#'
#'Using the default option, the GC normalization estimates the dependency between reads count and GC content only on the not methylated region,
#'This strategy can avoid the quantification biases for the modifications that biologically favor extream GC content.
#'
#'The GC normalization can results in an improvement of accuracy for most published m6A-seq data,
#'and it is particullarly recommended if you want to compare the data between different laboratory conditions.
#'
#'@details By default, the sequencing depth factor used is the column size factors obtained by \code{\link{estimateSeqDepth}},
#' if the column size factors are not provided, it will be estimated using the default method defined in \code{\link{estimateSeqDepth}}.
#'
#'@return a \code{summarizedExomePeak} object containing the feature specific size factors.
#'
#'@import SummarizedExperiment
#'@importFrom BSgenome getBSgenome
#'@docType methods
#'
#'@name normalizeGC
#'
#'@rdname normalizeGC
#'
#'@export

setMethod("normalizeGC",
          "SummarizedExomePeak",
                        function(sep,
                                 bsgenome = NULL,
                                 txdb = NULL,
                                 gene_anno_gff = NULL,
                                 fragment_length = 100,
                                 binding_length = 25,
                                 feature = c("background", "all"),
                                 qtnorm = TRUE,
                                 effective_GC = FALSE,
                                 drop_overlapped_genes = TRUE
                                 ) {
if(is.null(bsgenome)) {
stop("require BSgenome objects in GC size factor estimaton.")
}

if(!is.null(gene_anno_gff)) {
  txdb <- makeTxDbFromGFF(gene_anno_gff)
} else {
  if (is.null(txdb)) {
    stop("require transcript annotation in either GTF/GFF file or txdb object.")
  }

  if (!is(txdb, "TxDb")) {
    txdb <- makeTxDbFromUCSC(txdb)
  }
}

bsgenome <- getBSgenome(bsgenome)

#check if the sep object is abscent of the collumn wised size factors.
if(is.null(colData( sep )$sizeFactor)){
sep <- estimateSeqDepth(sep)
}

#CQN normalization with everything pooled
feature <- match.arg(feature)

#retrieve the vector of the exon level GC content from Genome.

elementMetadata( sep ) <- GC_content_over_grl(
                                      bsgenome = bsgenome,
                                      txdb = txdb,
                                      grl = rowRanges( sep ),
                                      fragment_length = fragment_length,
                                      binding_length = binding_length,
                                      drop_overlapped_genes = drop_overlapped_genes,
                                      effective_GC = effective_GC
                                 )

#Reserve any potential NA in the GC_content vector
GC_size_factors <- matrix(NA, nrow = nrow(sep), ncol = ncol(sep))

rownames(GC_size_factors) = rownames(sep)

GC_na_index <- is.na(elementMetadata( sep )$GC_content)


if(feature == "all"){
  Subindex = which( rowMeans(assay(sep)[!GC_na_index,]) > 50 )
} else {
  Subindex = which( rowMeans(assay(sep)[!GC_na_index,]) > 50 & grepl("control",rownames(sep))[!GC_na_index] )
}

if(!qtnorm) {

cqnObject <- suppressMessages( cqn(assay(sep)[!GC_na_index,],
                                     lengths = elementMetadata( sep )$feature_length[!GC_na_index],
                                     lengthMethod = "smooth",
                                     x = elementMetadata( sep )$GC_content[!GC_na_index],
                                     subindex = Subindex,
                                     sizeFactors = sep$sizeFactor,
                                     sqn = qtnorm,
                                     verbose = FALSE) )

GC_size_factors[!GC_na_index,] <- cqnObject$glm.offset

assays(sep)$GCsizeFactors <- GC_size_factors

} else {

cqnObject_IP <- suppressMessages( cqn(assay(sep)[!GC_na_index, sep$design_IP],
                                     lengths = elementMetadata( sep )$feature_length[!GC_na_index],
                                     lengthMethod = "smooth",
                                     x = elementMetadata( sep )$GC_content[!GC_na_index],
                                     subindex = Subindex,
                                     sizeFactors = sep$sizeFactor[sep$design_IP],
                                     sqn = qtnorm,
                                     verbose = FALSE) )

cqnObject_input <- suppressMessages( cqn(assay(sep)[!GC_na_index, !sep$design_IP],
                                      lengths = elementMetadata( sep )$feature_length[!GC_na_index],
                                      lengthMethod = "smooth",
                                      x = elementMetadata( sep )$GC_content[!GC_na_index],
                                      subindex = Subindex,
                                      sizeFactors = sep$sizeFactor[!sep$design_IP],
                                      sqn = qtnorm,
                                      verbose = FALSE) )

GC_size_factors[!GC_na_index,sep$design_IP] <- cqnObject_IP$glm.offset

GC_size_factors[!GC_na_index,!sep$design_IP] <- cqnObject_input$glm.offset

rm(cqnObject_IP, cqnObject_input)

assays(sep)$GCsizeFactors <- GC_size_factors

}

return(sep)

})

