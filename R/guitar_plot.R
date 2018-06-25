#' @title Plot the transcript distribution of the peaks
#'
#' @description  This function plot the transcript topology of RNA modification based on the provided SummarizedExomePeaks object.
#'
#' @details
#' If the SummarizedExomePeaks object contains quantification results for methylation, the significantly methylated peaks
#' with IP to input log2FC > 0 and DESeq2 Wald test padj < 0.05 will be plotted.
#'
#' If the SummarizedExomePeaks object contains quantification results for differential methylation, both the hyper methylation
#' and hypo methylation peaks with DESeq2 Wald test p values < 0.05 will be plotted.

#' @param sep a SummarizedExomePeaks object.
#' @param txdb a txdb object containing the transcript annotation.
#' @param save_pdf_prefix if provided, a pdf file with the given name will be saved under the current directory.
#' @param include_control_regions a logical indicating whether to include the control regions or not.
#' @param guitar_coordinate optional, the guitar coordinate of the transcript annotation.
#'
#' @return a ggplot object
#'
#' @importFrom Guitar GuitarPlot makeGuitarCoordsFromTxDb
#' @export
guitar_plot <- function(sep,
                        txdb = NULL,
                        save_pdf_prefix = NULL,
                        include_control_regions = TRUE,
                        guitar_coordinate = NULL) {

if(sum(grepl("meth",rownames(sep$SE))) < 10 ) {
  stop("guitar plot function cannot be performed for total peaks number < 10.")
}

stopifnot(!(is.null(txdb) & is.null(guitar_coordinate)))

if(is.null(guitar_coordinate)){
  guitar_coordinate = makeGuitarCoordsFromTxDb(txdb)
}

  #first check whether the input object contains quantification result.
    #if so, we need only plot the meth peaks and control peaks.
if(is.null(sep$DESeq2Result)){
row_grl <- rowRanges( sep$SE )

gr_list <- list(
  peaks = row_grl[grepl("meth",names(row_grl) )],
  control = row_grl[grepl("control",names(row_grl) )]
)

if(!include_control_regions){
  gr_list <- gr_list[-2]
}

suppressWarnings(
              GuitarPlot(
                 gfeatures = gr_list,
                 GuitarCoordsFromTxDb = guitar_coordinate,
                 saveToPDFprefix = ifelse(is.null(save_pdf_prefix),NA,save_pdf_prefix)
                        )
                )

} else {
  #In case of the collumn design contains only methylation
  if(!any(sep$SE$design_Treatment)){

   final_alpha  <- decision_deseq2(Inf_RES = sep$DESeq2Result,
                                   Padj_cut = 0.05,
                                   log2FC_cut = 0,
                                    Min_mod = 1000)$Cut_Val_expected

   if(is.na(final_alpha)) {final_alpha = 1}

   row_grl <- rowRanges( sep$SE )

   gr_list <- list(meth_peaks = row_grl[grepl("meth",rownames(sep$SE))][sep$DESeq2Result$padj < final_alpha],
        control = row_grl[grepl("control",names(row_grl) )]
   )

   if(!include_control_regions){
     gr_list <- gr_list[-2]
   }

   suppressWarnings(

                  Guitar::GuitarPlot(
                  gfeatures = gr_list,
                  GuitarCoordsFromTxDb = guitar_coordinate,
                  saveToPDFprefix = ifelse(is.null(save_pdf_prefix),NA,save_pdf_prefix)
                )

            )


  } else {
    #Finally, the guitar plot will be generated for differential methylation results.

    decision_table  <- decision_deseq2(Inf_RES = sep$DESeq2Result,
                                    P_cut = 0.05,
                                    log2FC_cut = 0,
                                    Min_mod = 1000,
                                    Exp_dir = "hyper")

    if(is.na(decision_table$Cut_Val_ctrl)) decision_table$Cut_Val_ctrl = 1

    if(is.na(decision_table$Cut_Val_expected)) decision_table$Cut_Val_expected = 1

    row_grl <- rowRanges( sep$SE )

    indx_hyper <- which( sep$DESeq2Result$pvalue < decision_table$Cut_Val_expected &
                         sep$DESeq2Result$log2FoldChange > 0)

    indx_hypo <- which( sep$DESeq2Result$pvalue < decision_table$Cut_Val_ctrl &
                           sep$DESeq2Result$log2FoldChange < 0)

    gr_list <- list(hyperMeth = row_grl[grepl("meth",rownames(sep$SE))][indx_hyper],
                    hypoMeth = row_grl[grepl("meth",rownames(sep$SE))][indx_hypo]
    )

    suppressWarnings(

      GuitarPlot(
        gfeatures = gr_list,
        GuitarCoordsFromTxDb = guitar_coordinate,
        saveToPDFprefix = ifelse(is.null(save_pdf_prefix),NA,save_pdf_prefix)
      )

    )

}

}

}
