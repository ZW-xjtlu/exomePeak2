#' @title Calculate the decision table for a DESeq2 result.
#'
#' @description \code{Decision_dsresult} is an internal function used to summary the cut-off and the number of positive results used for DESeq2 result..
#'
#' @param Inf_RES A \code{DESeqResults} or similar object that contains the result statistics for either methylation or differential methylation.
#' @param log2FC_cut The log2 fold change cutoff of the inference result, default setting is 0.
#' @param P_cut A numeric value between 0 to 1, indicating the p value cut off of the Wald test defined by DESeq2 (or defined by QNB), it will be neglected if \code{Padj_cut} is not NULL.
#' @param Padj_cut A numeric value between 0 to 1, indicating the fdr cut off of the Wald test defined by DESeq2 (or defined by QNB).
#' @param Min_mod Minimum number of features returned, when this is smaller than the cut-off results, additional features are called by the order of p values.
#' @param Exp_dir This parameter is filled when making decisions on differential methylation, it could be "hyper", "hypo", and "both".
#'
#' @return A \code{data.frame} object indicating the column and cut-off value used for desicion, also it includes the number of positive sites in both directions based on the decision.
#'
#' @export
#'
decision_deseq2 <- function(Inf_RES,
                            log2FC_cut = 0,
                            P_cut = 0.05,
                            Padj_cut = NULL,
                            Min_mod = 1000,
                            Exp_dir = c("hyper","hypo","both")){

  Exp_dir <- match.arg(Exp_dir)

  Inf_RES <- Inf_RES[!(is.na(Inf_RES$log2FoldChange) | is.na(Inf_RES$pvalue)),]

  Inf_RES <- Inf_RES[abs(Inf_RES$log2FoldChange) > log2FC_cut,]

  result_df <- data.frame(
    log2FC_cut = log2FC_cut,
    Cut_By_ctrl = "pvalue" ,
    Cut_By_expected = "pvalue" ,
    Cut_Val_ctrl = 0 ,
    Cut_Val_expected = 0 ,
    Discoveries_ctrl = 0 ,
    Discoveries_expected = 0,
    Expected_dir = "> 0"
  )

  if(!is.null(Exp_dir)) {
    if(Exp_dir=="hypo") {result_df$Expected_dir = "< 0"
     } else {
    if(Exp_dir=="both") result_df$Expected_dir = "!= 0"
  }
  }

  #Define index of the expected direction.
  EXP_idx =  eval(parse(text = paste0( "Inf_RES$log2FoldChange ",result_df$Expected_dir) ))

  #First let's do regular cut-off

  if( !is.null( Padj_cut ) ) {
    Inf_RES$padj[is.na(Inf_RES$padj)] = 1
    result_df$Cut_By_ctrl = "padj"
    result_df$Cut_By_expected = "padj"
    result_df$Cut_Val_ctrl = Padj_cut
    result_df$Cut_Val_expected = Padj_cut
    result_df$Discoveries_ctrl = sum(!EXP_idx & Inf_RES$padj < Padj_cut)
    result_df$Discoveries_expected = sum(EXP_idx & Inf_RES$padj < Padj_cut)
  } else {
    result_df$Cut_Val_ctrl = P_cut
    result_df$Cut_Val_expected = P_cut
    result_df$Discoveries_ctrl = sum(!EXP_idx & Inf_RES$pvalue < P_cut)
    result_df$Discoveries_expected = sum(EXP_idx & Inf_RES$pvalue < P_cut)
  }

  #Second, we will consider the case while the minimum number is not met.

  if(Exp_dir == "both") {

    if(result_df$Discoveries_expected < Min_mod) {
      result_df$Cut_By_expected = "pvalue"

      if( nrow(Inf_RES) <= Min_mod ) {
        result_df$Cut_Val_expected = sort( Inf_RES$pvalue )[ nrow(Inf_RES) + 1 ]
        result_df$Discoveries_expected = sum(Inf_RES$pvalue < result_df$Cut_Val_expected)
      } else {
        result_df$Cut_Val_expected = sort( Inf_RES$pvalue )[ Min_mod + 1 ]
        result_df$Discoveries_expected = sum(Inf_RES$pvalue < result_df$Cut_Val_expected)
      }
    }

  } else {
    if(result_df$Discoveries_ctrl < Min_mod) {

      result_df$Cut_By_ctrl = "pvalue"

      if( sum(!EXP_idx) <= Min_mod ) {
        result_df$Cut_Val_ctrl = sort( Inf_RES$pvalue[!EXP_idx])[ sum(!EXP_idx) + 1 ]
        result_df$Discoveries_ctrl = sum(Inf_RES$pvalue[!EXP_idx] < result_df$Cut_Val_ctrl)
      } else {
        result_df$Cut_Val_ctrl = sort( Inf_RES$pvalue[!EXP_idx])[ Min_mod + 1 ]
        result_df$Discoveries_ctrl = sum(Inf_RES$pvalue[!EXP_idx] < result_df$Cut_Val_ctrl)
      }
    }

    if(result_df$Discoveries_expected < Min_mod) {
      result_df$Cut_By_expected = "pvalue"
      if( sum(EXP_idx) <= Min_mod ) {
        result_df$Cut_Val_expected = sort( Inf_RES$pvalue[EXP_idx])[ sum(EXP_idx) + 1 ]
        result_df$Discoveries_expected = sum(Inf_RES$pvalue[EXP_idx] < result_df$Cut_Val_expected)
      } else {
        result_df$Cut_Val_expected = sort( Inf_RES$pvalue[EXP_idx])[ Min_mod + 1 ]
        result_df$Discoveries_expected = sum(Inf_RES$pvalue[EXP_idx] < result_df$Cut_Val_expected)
      }
    }
  }

  return(result_df)
}
