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
#' @docType methods
#'
#' @name plotGuitar
#'
#' @rdname plotGuitar
#' @export

setMethod("plotGuitar",
          "SummarizedExomePeak",
               function(sep,
                        txdb = NULL,
                        save_pdf_prefix = NULL,
                        include_control_regions = TRUE,
                        guitar_coordinate = NULL) {

if(sum(grepl("meth",rownames(sep))) < 10 ) {
  stop("guitar plot function cannot be performed for total peaks number < 10.")
}

stopifnot(!(is.null(txdb) & is.null(guitar_coordinate)))

if(is.null(guitar_coordinate)){
  guitar_coordinate = makeGuitarCoordsFromTxDb(txdb)
}

  #first check whether the input object contains quantification result.
    #if so, we need only plot the meth peaks and control peaks.
if(nrow(DESeq2Results(sep)) == 0){
row_grl <- rowRanges( sep )

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
  if(!any(sep$design_Treatment)){

    row_grl <- rowRanges( sep )

    indx_sig <- which( DESeq2Results(sep)$padj < .05 & DESeq2Results(sep)$log2FoldChange > 0 )

    gr_lab = "mod padj < .05"

    if( length(indx_sig) < floor( sum(grepl("meth_", rownames(sep))) * 0.01 ) ){

      indx_sig <- which( DESeq2Results(sep)$pvalue < .05 & DESeq2Results(sep)$log2FoldChange > 0 )

      gr_lab = "mod p < .05"

    }

    gr_list <- list(meth_peaks = row_grl[ grepl("meth", rownames(sep)) ][indx_sig],
                    control = row_grl[ grepl("control", names(row_grl)) ]
    )

    names(gr_list)[1] <- gr_lab

    rm(indx_sig, gr_lab)

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

    row_grl <- rowRanges( sep )

    indx_hyper <- which( DESeq2Results(sep)$padj < .05 &
                           DESeq2Results(sep)$log2FoldChange > 0)

    indx_hypo <- which( DESeq2Results(sep)$padj < .05 &
                          DESeq2Results(sep)$log2FoldChange < 0)

    list_names <- c("hyper padj < .05", "hypo padj < .05")

    min_positive <- floor(sum(grepl("meth_", rownames(sep))) * 0.005)

    if(length(indx_hyper) + length(indx_hypo) < min_positive){

      indx_hyper <- which( DESeq2Results(sep)$padj < .05 &
                             DESeq2Results(sep)$log2FoldChange > 0)

      indx_hypo <- which( DESeq2Results(sep)$padj < .05 &
                            DESeq2Results(sep)$log2FoldChange < 0)

      list_names <- c("hyper p < .05", "hypo p < .05")

    }

    gr_list <- list(hyperMeth = row_grl[grepl("meth",rownames(sep))][indx_hyper],
                    hypoMeth = row_grl[grepl("meth",rownames(sep))][indx_hypo]
    )

    names(gr_list) <- list_names

    suppressWarnings(

      GuitarPlot(
        gfeatures = gr_list,
        GuitarCoordsFromTxDb = guitar_coordinate,
        saveToPDFprefix = ifelse(is.null(save_pdf_prefix),NA,save_pdf_prefix)
      )

    )

}

}

})
