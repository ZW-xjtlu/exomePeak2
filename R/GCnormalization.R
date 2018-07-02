#'@title Estimation of the GC content normalization factors.
#'@param sep a \code{summarizedExomePeak} object.
#'
#'@param bsgenome a \code{\link{BSgenome}} object for the genome sequence, or it could be the name of this reference genome specified in a way that is accepted by the getBSgenome function defined in the BSgenome software package
#'
#'@param feature the features used in the GC effect estimation, can be "background" and "all".
#'If "all" is choosed, the GC effect function will be estimated using both the methylated and the background region,
#'this choice will force the resulting modification signals independent of GC content.
#'However,it is possibably biased if the modification is biologically dependent on GC content;
#'Default "background".
#'
#'@param qtnorm a logical indicating wheather perform quantile normalization after the GC effect estimation.
#'Quantile normalization can correct systematic effect in RNA-seq count data.
#'However, it may result in errors due to the biological difference between the distributions of IP and input.
#'
#'@param fragment_length the expected fragment length of the sequencing library.
#'The widths of the row features for quantifying GC content will be re-sized into the fragment length if it falls bellow it.
#'
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
#'@details By default, the sequencing depth factor used is the column size factors obtained by \code{\link{estimate_size_factors}},
#' if the column size factors are not provided, it will be estimated using the default method defined in \code{\link{estimate_size_factors}}.
#'
#'@return a \code{summarizedExomePeak} object containing the feature specific size factors.
#'
#'@import SummarizedExperiment
#'@importFrom BSgenome getBSgenome
#'@docType methods
#'
#'@name GCnormalization
#'
#'@rdname GCnormalization
#'
#'@export

setMethod("GCnormalization",
          "SummarizedExomePeak",
                        function(sep,
                                 bsgenome = "hg19",
                                 feature = c("background","all"),
                                 qtnorm = FALSE,
                                 fragment_length = 100
                                 ) {

stopifnot(!(is.null(bsgenome)))

if(is.character(bsgenome)) {
  bsgenome <- getBSgenome(bsgenome)
}

#check if the sep object is abscent of the collumn wised size factors.
if(is.null(colData( sep )$sizeFactor)){
sep <- estimateSeqDepth(sep)
}

#CQN normalization with everything pooled
feature <- match.arg(feature)

#retrieve the vector of GC content from Genome.
CQN_parameter <- GC_content_over_grl(bsgenome = bsgenome,
                                  grl = rowRanges( sep ),
                                  fragment_length = fragment_length)

#Reserve any potential NA in the GC_content vector
GC_size_factors <- matrix(NA,nrow = nrow(sep),ncol = ncol(sep))
rownames(GC_size_factors) = rownames(sep)
GC_na_index <- is.na(CQN_parameter$GC_content)


if(feature == "all"){
  Subindex = which( rowMeans(assay(sep)) > 50 )
} else {
  Subindex = which( rowMeans(assay(sep)) > 50 & grepl("control",rownames(sep)) )
}

cqnObject <- suppressMessages( cqn(assay(sep)[!GC_na_index,],
                                   lengths = CQN_parameter$Indx_length[!GC_na_index],
                                   lengthMethod = "smooth",
                                   x = CQN_parameter$GC_content[!GC_na_index],
                                   subindex = Subindex,
                                   sizeFactors = sep$sizeFactor,
                                   sqn = qtnorm,
                                   verbose = FALSE) )

cqnOffset <- cqnObject$glm.offset
normFactors <- exp(cqnOffset)
normFactors <- normFactors / exp(rowMeans(log(normFactors)))

GC_size_factors[!GC_na_index,] <- normFactors

GCsizeFactors( sep )<- GC_size_factors

return(sep)

})
