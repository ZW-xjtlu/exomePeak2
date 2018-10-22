#' @title Calculate the decision table for a DESeq result.
#'
#' @description \code{Decision_dsresult} is an internal function used to summary the cut-off and the number of positive results used for DESeq result..
#'
#' @param res A \code{DESeqResults} or similar object that contains the result statistics for either methylation or differential methylation.
#' @param log2FC_cut The log2 fold change cutoff of the inference result, default setting is 0.
#' @param p_cut A numeric value between 0 to 1, indicating the p value cut off of the Wald test defined by DESeq, it will be neglected if \code{padj_cut} is not NULL.
#' @param padj_cut A numeric value between 0 to 1, indicating the fdr cut off of the Wald test defined by DESeq.
#' @param min_mod Minimum number of features returned, when this is smaller than the cut-off results, additional features are called by the order of p values.
#'
#' @return A \code{data.frame} object indicating the column and cut-off value used for desicion, also it includes the number of positive sites in both directions based on the decision.
#'
#' @export
decision_deseq <- function(res,
                           log2FC_cut,
                           p_cut,
                           padj_cut,
                           min_mod) {

  res <- res[!(is.na(res$log2FoldChange) | is.na(res$pval)),]

  res <- res[abs(res$log2FoldChange) > log2FC_cut,]

  result_df <- data.frame(
    log2FC_cut = log2FC_cut,
    Cut_By_ctrl = "pval" ,
    Cut_By_expected = "pval" ,
    Cut_Val_ctrl = 0 ,
    Cut_Val_expected = 0 ,
    Discoveries_ctrl = 0 ,
    Discoveries_expected = 0,
    Expected_dir = "> 0"
  )

  #Define index of the expected direction.
  EXP_idx =  eval(parse(text = paste0( "res$log2FoldChange ",result_df$Expected_dir) ))

  #First let's do regular cut-off
  if( !is.null( padj_cut ) ) {
    res$padj[is.na(res$padj)] = 1
    result_df$Cut_By_ctrl = "padj"
    result_df$Cut_By_expected = "padj"
    result_df$Cut_Val_ctrl = padj_cut
    result_df$Cut_Val_expected = padj_cut
    result_df$Discoveries_ctrl = sum(!EXP_idx & res$padj < padj_cut)
    result_df$Discoveries_expected = sum(EXP_idx & res$padj < padj_cut)
  } else {
    result_df$Cut_Val_ctrl = P_cut
    result_df$Cut_Val_expected = P_cut
    result_df$Discoveries_ctrl = sum(!EXP_idx & res$pval < P_cut)
    result_df$Discoveries_expected = sum(EXP_idx & res$pval < P_cut)
  }

  #Second, we will consider the case while the minimum number is not met.
  if(result_df$Discoveries_ctrl < min_mod) {
    result_df$Cut_By_ctrl = "pval"
    result_df$Cut_Val_ctrl = sort( res$pval[!EXP_idx])[min_mod+1]
    result_df$Discoveries_ctrl = sum(res$pval[!EXP_idx] < result_df$Cut_Val_ctrl)
  }

  if(result_df$Discoveries_expected < min_mod) {
    result_df$Cut_By_expected = "pval"
    result_df$Cut_Val_expected = sort( res$pval[EXP_idx])[min_mod+1]
    result_df$Discoveries_expected = sum(res$pval[EXP_idx] < result_df$Cut_Val_expected)
  }

  return(result_df)
}
