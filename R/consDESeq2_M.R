#' @title Calculate the index for across sample modification consistency in a DESeqDataSet object
#'
#' @param dds a DESeqDataSet object.
#'
#' @param consistent_log2FC_cutoff a \code{numeric} for the modification log2 fold changes cutoff in the peak consisency calculation; default = 1.
#'
#' @param consistent_fdr_cutoff a \code{numeric} for the BH adjusted C-test p values cutoff in the peak consistency calculation; default { = 0.05}. Check \code{\link{ctest}}.
#'
#' @param alpha a \code{numeric} for the binomial quantile used in the consitent peak filter; default\code{ = 0.05}.
#'
#' @param p0 a \code{numeric} for the binomial proportion parameter used in the consistent peak filter; default \code{= 0.8}.
#'
#' For a peak to be consistently methylated, the minimum number of significant enriched replicate pairs is defined as the 1 - alpha quantile of a binomial distribution with p = p0 and N = number of possible pairs between replicates.
#'
#' The consistency defined in this way is equivalent to the rejection of an exact binomial test with null hypothesis of p < p0 and N = replicates number of IP * replicates number of input.
#'
#' @details The minimum consistent number cutcoff is defined by 1-alpha quantile of a binomial distribution with probability of success = p, and number of trials = number of possible pairs between replicates.
#' This is equivalent to the rejection of an exact binomial test with null hypothesis of p < 0.8 and N = number of possible pairs.
#' @import SummarizedExperiment
#'
#'@return a logical index for the consistently modified rows in the DESeqDataSet.
#'
consDESeq2_M <- function(dds,
                        consistent_log2FC_cutoff = 1,
                        consistent_fdr_cutoff = 0.05,
                        p0 = 0.8,
                        alpha = 0.05){

#normalize the count with size factors
if(is.null(normalizationFactors(dds))){
  sizeFactor_M <- matrix(dds$sizeFactor,nrow = nrow(dds), ncol = ncol(dds), byrow = T)
}else{
  sizeFactor_M <- normalizationFactors(dds)
}

N = sum(!dds$design_Treatment & dds$design_IP == "IP") * sum(!dds$design_Treatment & dds$design_IP == "input")

cut = qbinom(1-alpha, size = N, prob = p0)

consist_count = rep(0, nrow(dds))

for (i in which(!dds$design_Treatment & dds$design_IP == "IP")){
  for( j in which(!dds$design_Treatment & dds$design_IP == "input")){
     test_result <- ctest(assay(dds)[,i],
                          assay(dds)[,j],
                          sizeFactor_M[,i],
                          sizeFactor_M[,j],
                          fold = 1)
     consist_count <- consist_count + as.numeric(test_result$fdr < consistent_fdr_cutoff & test_result$log2FC > consistent_log2FC_cutoff)
  }
}

consist_indx <- consist_count >= cut

rm(consist_count,N,cut)

if(any(dds$design_Treatment)){

  N = sum(dds$design_Treatment & dds$design_IP == "IP") * sum(dds$design_Treatment & dds$design_IP == "input")

  cut = qbinom(1-alpha,size = N, prob = p0)

  consist_countT = rep(0, nrow(dds))

  for (i in which(dds$design_Treatment & dds$design_IP == "IP")){
    for( j in which(dds$design_Treatment & dds$design_IP == "input")){

      test_result <- ctest(assay(dds)[,i],
                           assay(dds)[,j],
                           sizeFactor_M[,i],
                           sizeFactor_M[,j],
                           fold = 1,
                           alpha = consistent_fdr_cutoff)

      consist_countT <- consist_countT + as.numeric(test_result$fdr < consistent_fdr_cutoff & test_result$log2FC > consistent_log2FC_cutoff)
    }
  }
consist_indx = consist_indx | consist_countT >= cut
}

return(consist_indx)
}
