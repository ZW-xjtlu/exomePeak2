#'@title Estimate Normalization Factors for GC content Bias Correction.
#'
#'@param sep a \code{\link{summarizedExomePeak}} object returned by \code{\link{exomePeak2}} or \code{\link{exomePeakCalling}}.
#'
#'@param txdb a \code{\link{TxDb}} object for the transcript annotation,
#' If the \code{TxDb} object is not available, it could be a \code{character} string of the UCSC genome name which is acceptable by \code{\link{makeTxDbFromUCSC}}, example: \code{"hg19"}.
#'
#'@param bsgenome a \code{\link{BSgenome}} object for the genome reference,
#' If the \code{BSgenome} object is not available, it could be a \code{character} string of the UCSC genome name which is acceptable by \code{\link{getBSgenome}}, example: \code{"hg19"}.
#'
#'@param fragment_length a positive integer number for the expected fragment length in nucleotides; default \code{= 100}.
#'
#'@param binding_length a positive integer number for the expected binding length of the anti-modification antibody in IP samples; default \code{= 25}.
#'
#'@param feature a \code{character} specifies the region used in the GC effect estimation, can be one in \code{c("all", "background")}; default \code{"all"}.
#'
#'\describe{
#'
#'  \item{\strong{\code{all}}}{
#'  The GC content and the parameters for quantile normalization will be estimated on all regions.
#'  }
#'
#'  \item{\strong{\code{background}}}{
#'  The GC content linear effect and the parameters for quantile normalization will be estimated on the background regions.
#'  }
#'
#' }
#'
#'@param qtnorm a \code{logical} of whether to perform subset quantile normalization after the GC content linear effect correction； default \code{= TRUE}.
#'
#' If \code{qtnorm = TRUE}, subset quantile normalization will be applied within the IP and input samples seperately to account for the inherent differences between the marginal distributions of IP and input samples.
#'
#'@param effective_gc a \code{logical} of whether to calculate the effective GC content weighted by the fragment alignment probabilities; default \code{= FALSE}.
#'
#'@description \code{normalizeGC} estimates the feature specific size factors in order to reduce the technical variation during modification peak statistics quantification.
#'
#'
#'@details
#'PCR amplication bias related to GC content is a major source of technical variation in RNA-seq.
#'The GC content biases are usually correlated within the same laboratory environment, and this will result in the batch effect between different studies.
#'
#'The GC content normalization can result in an improvement of peak accuracy for most published m6A-seq data,
#'and it is particullarly recommended if you want to compare the quantifications on methylation levels between different laboratory conditions.
#'
#'@return a \code{summarizedExomePeak} object with the updated slot \code{GCsizeFactors}.
#'
#'@examples
#'
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#'
#'aln <- scanMeripBAM(
#'bam_ip = c("IP_rep1.bam",
#'           "IP_rep2.bam",
#'           "IP_rep3.bam"),
#'bam_input = c("input_rep1.bam",
#'              "input_rep2.bam",
#'              "input_rep3.bam"),
#'paired_end = TRUE
#')
#'
#'sep <- exomePeakCalling(merip_bams = aln,
#'                        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                        bsgenome = Hsapiens)
#'
#'sep <- normalizeGC(sep)
#'
#'@import SummarizedExperiment
#'@import cqn
#'@importFrom BSgenome getBSgenome
#'@docType methods
#'@seealso \code{\link{estimateSeqDepth}}
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
                                 gff_dir = NULL,
                                 fragment_length = 100,
                                 binding_length = 25,
                                 feature = c("all","background"),
                                 qtnorm = TRUE,
                                 effective_GC = FALSE
                                 ) {
feature <- match.arg(feature)

#check if the sep object is abscent of the collumn wised size factors.
if(is.null(colData( sep )$sizeFactor)){
  sep <- estimateSeqDepth(sep)
}

if(any(is.null(elementMetadata( sep )$GC_content),
          is.null(elementMetadata( sep )$feature_length))) {

if(is.null(bsgenome)) {
stop("Require BSgenome objects in GC size factor estimaton.")
}

if(!is.null(gff_dir)) {
  txdb <- makeTxDbFromGFF(gff_dir)
} else {
  if (is.null(txdb)) {
    stop("require transcript annotation in either GTF/GFF file or txdb object.")
  }

  if (!is(txdb, "TxDb")) {
    txdb <- makeTxDbFromUCSC(txdb)
  }
}

bsgenome <- getBSgenome(bsgenome)

#retrieve the vector of the exon level GC content from Genome.

elementMetadata( sep ) <- GC_content_over_grl(
                                      bsgenome = bsgenome,
                                      txdb = txdb,
                                      grl = rowRanges( sep ),
                                      fragment_length = fragment_length,
                                      binding_length = binding_length,
                                      effective_GC = effective_GC
                                 )

}

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

cqnObject <- quiet( suppressMessages(
                    suppressWarnings(
                                     cqn(assay(sep)[!GC_na_index,],
                                     lengths = elementMetadata( sep )$feature_length[!GC_na_index],
                                     lengthMethod = "fixed",
                                     x = elementMetadata( sep )$GC_content[!GC_na_index],
                                     subindex = Subindex,
                                     sizeFactors = sep$sizeFactor,
                                     sqn = qtnorm,
                                     verbose = FALSE) )
                    )
                    )

GC_size_factors[!GC_na_index,] <- cqnObject$glm.offset

assays(sep)$GCsizeFactors <- GC_size_factors

} else {

cqnObject_IP <- quiet(  suppressMessages(
                        suppressWarnings(
                                     cqn(assay(sep)[!GC_na_index, sep$design_IP],
                                     lengths = elementMetadata( sep )$feature_length[!GC_na_index],
                                     lengthMethod = "fixed",
                                     x = elementMetadata( sep )$GC_content[!GC_na_index],
                                     subindex = Subindex,
                                     sizeFactors = sep$sizeFactor[sep$design_IP],
                                     sqn = qtnorm,
                                     verbose = FALSE)
                                     )
                               )
                        )

cqnObject_input <-  quiet( suppressMessages(
                           suppressWarnings(
                                      cqn(assay(sep)[!GC_na_index, !sep$design_IP],
                                      lengths = elementMetadata( sep )$feature_length[!GC_na_index],
                                      lengthMethod = "fixed",
                                      x = elementMetadata( sep )$GC_content[!GC_na_index],
                                      subindex = Subindex,
                                      sizeFactors = sep$sizeFactor[!sep$design_IP],
                                      sqn = qtnorm,
                                      verbose = FALSE)
                                      )
                                )
                           )

GC_size_factors[!GC_na_index,sep$design_IP] <- cqnObject_IP$glm.offset

GC_size_factors[!GC_na_index,!sep$design_IP] <- cqnObject_input$glm.offset

rm(cqnObject_IP, cqnObject_input)

assays(sep)$GCsizeFactors <- GC_size_factors

}

return(sep)

})

