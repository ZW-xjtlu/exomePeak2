#'@title Calculate the index for across sample differential modification consistency in a DESeqDataSet object
#'@param dds a DESeqDataSet object.
#'@param consistent_log2FC_cutoff a numeric for the cutoff of differential modification log2 fold changes; default = 1.
#'@param p the binomial parameter p used by the binomial test; default = 0.7.
#'@param alpha the significant level of the binomial test; default = 0.05.
#'@details The minimum consistent number cutcoff is defined by 1-alpha quantile of a binomial distribution with probability of success = p, and number of trials = number of possible pairs between replicates.
#'This is equivalent to the rejection of an exact binomial test with null hypothesis of p < 0.8 and N = number of possible pairs.
#'@import SummarizedExperiment
#'
#'@return a logical index for the consistently modified rows in the DESeqDataSet.
#'@keywords internal
#'
consDESeq2_DM <- function(dds,
                         consistent_log2FC_cutoff = 1,
                         p = 0.7,
                         alpha = 0.05){

  stopifnot(any(dds$design_Treatment))

  #normalize the count with size factors
  if(is.null(normalizationFactors(dds))){
    abundence <- assay(dds)/dds$sizeFactor
  }else{
    abundence <- assay(dds)/normalizationFactors(dds)
  }

  N = sum(!dds$design_Treatment & dds$design_IP == "IP") *
    sum(!dds$design_Treatment & dds$design_IP == "input") *
    sum(dds$design_Treatment & dds$design_IP == "IP") *
    sum(dds$design_Treatment & dds$design_IP == "input")

  cut = qbinom(1-alpha,size = N, prob = p)

  consist_count_up = rep(0, nrow(dds))

  consist_count_down = rep(0, nrow(dds))

  # Both the consistent up regulation and the consistent down regulation will be regarded as consistent differential modification

  for (i in which(!dds$design_Treatment & dds$design_IP == "IP")){
    for( j in which(!dds$design_Treatment & dds$design_IP == "input")){
      for ( k in which(dds$design_Treatment & dds$design_IP == "IP")){
        for ( l in which(dds$design_Treatment & dds$design_IP == "input")){
          diff_LFC = (abundence[,k]+0.001/abundence[,l]+0.001)/(abundence[,i]+0.001/abundence[,j]+0.001)
          consist_count_up = consist_count_up + as.numeric(diff_LFC > exp(consistent_log2FC_cutoff))
          consist_count_down = consist_count_down + as.numeric(diff_LFC < exp(consistent_log2FC_cutoff))
        }
      }
    }
  }

  consist_indx <- consist_count_up >= cut | consist_count_down >= cut

  rm(consist_count_up,consist_count_down,diff_LFC,N,cut)

  return(consist_indx)
}
