#' @title Plot the exon length distribution of the peaks
#'
#' @description This function plot the distribution of the exon length for peaks containing exons.
#' @details
#' If the SummarizedExomePeaks object contains quantification results for methylation, the significantly methylated peaks
#' with IP to input log2FC > 0 and DESeq2 Wald test padj < 0.05 will be plotted .
#'
#' If the SummarizedExomePeaks object contains quantification results for differential methylation, both the hyper methylation
#' and hypo methylation peaks with DESeq2 Wald test p values < 0.05 will be plotted.
#'
#' @param sep a SummarizedExomePeaks object.
#' @param txdb a txdb object containing the transcript annotation.
#' @param save_pdf_prefix if provided, a pdf file with the given name will be saved under the current working directory.
#' @param include_control_regions a logical indicating whether to include the control regions or not.
#'
#' @return a ggplot object
#'
#' @import Guitar
#' @import SummarizedExperiment
#'
#' @docType methods
#'
#' @name plotExonLength
#'
#' @rdname plotExonLength
#'
#' @export
setMethod("plotExonLength",
          "SummarizedExomePeak",
                  function(sep,
                           txdb = NULL,
                           save_pdf_prefix = NULL,
                           include_control_regions = TRUE) {

if(sum(grepl("meth",rownames(sep))) < 10 ) {
    stop("exon length plot cannot be performed for total peaks number < 10.")
}

  stopifnot(!is.null(txdb))

  #first check whether the input object contains any quantification result.

  #if so, we need only plot the meth peaks and control peaks.

  if(is.null(DESeq2Results(sep))){
    row_grl <- rowRanges( sep )

    gr_list <- list(
      peaks = row_grl[grepl("meth",names(row_grl) )],
      control = row_grl[grepl("control",names(row_grl) )]
    )

    if(!include_control_regions){
      gr_list <- gr_list[-2]
    }

    suppressWarnings(

      ExonPlot(
        gfeatures = gr_list,
        txdb = txdb,
        save_pdf_prefix = save_pdf_prefix
      )

    )

  } else {
    #In case of the collumn design contains only methylation
    if(!any(sep$design_Treatment)){

      final_alpha  <- decision_deseq2(Inf_RES = DESeq2Results(sep),
                                      Padj_cut = 0.05,
                                      log2FC_cut = 0,
                                      Min_mod = 1000)$Cut_Val_expected

      if(is.na(final_alpha)) {final_alpha = 1}

      row_grl <- rowRanges( sep )

      gr_list <- list(meth_peaks = row_grl[grepl("meth",rownames(sep))][DESeq2Results(sep)$padj < final_alpha],
                      control = row_grl[grepl("control",names(row_grl) )]
      )

      if(!include_control_regions){
        gr_list <- gr_list[-2]
      }

      suppressWarnings(

        ExonPlot(
          gfeatures = gr_list,
          txdb = txdb,
          save_pdf_prefix = save_pdf_prefix
        )

      )


    } else {
      #Finally, the guitar plot will be generated for differential methylation results.

      decision_table  <- decision_deseq2(Inf_RES = DESeq2Results(sep),
                                         P_cut = 0.05,
                                         log2FC_cut = 0,
                                         Min_mod = 1000,
                                         Exp_dir = "hyper")

      if(is.na(decision_table$Cut_Val_ctrl)) decision_table$Cut_Val_ctrl = 1

      if(is.na(decision_table$Cut_Val_expected)) decision_table$Cut_Val_expected = 1

      row_grl <- rowRanges( sep )

      indx_hyper <- which( DESeq2Results(sep)$pvalue < decision_table$Cut_Val_expected &
                             DESeq2Results(sep)$log2FoldChange > 0)

      indx_hypo <- which( DESeq2Results(sep)$pvalue < decision_table$Cut_Val_ctrl &
                            DESeq2Results(sep)$log2FoldChange < 0)

      gr_list <- list(hyperMeth = row_grl[grepl("meth",rownames(sep))][indx_hyper],
                      hypoMeth = row_grl[grepl("meth",rownames(sep))][indx_hypo]
      )

      suppressWarnings(

        ExonPlot(
          gfeatures = gr_list,
          txdb = txdb,
          save_pdf_prefix = save_pdf_prefix
        )

      )

  }

}
})
