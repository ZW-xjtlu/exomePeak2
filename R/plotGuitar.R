#' @title Plot the distribution of the peaks/sites on travis coordinate.
#'
#' @description \code{plotGuitar} provides visualization of the peaks/sites' distribution on travis coordinate for the \code{\link{SummarizedExomePeaks}} object.
#'
#' @details
#' If the SummarizedExomePeaks object contains quantification results for modification, the significantly modified peaks
#' with IP to input log2FC > 0 and DESeq2 Wald test padj < 0.05 will be plotted.
#'
#' If the SummarizedExomePeaks object contains quantification results for differential modification, both the hyper modification
#' and hypo modification peaks with DESeq2 Wald test p values < 0.05 will be plotted.

#' @param sep a \code{\link{SummarizedExomePeaks}} object.
#' @param txdb a \code{\link{TxDb}} object containing the transcript annotation.
#' @param save_pdf_prefix a \code{character}, if provided, a pdf file with the given name will be saved under the current directory.
#' @param include_control_regions a \code{logical} for whether to plot the control regions together with the modification regions.
#' @param guitar_coordinate optional, the guitar coordinate of the transcript annotation.
#' @param save_dir optional, a \code{character} indicating the directory to save the plot; default ".".
#'
#' @return a \code{ggplot} object
#'
#' @examples
#'\dontrun{
#'sep #a SummarizedExomePeak object
#'txdb #a TxDb object
#'plotExonLength(sep,txdb)
#'}
#'
#' @docType methods
#'
#' @name plotGuitar
#'
#' @rdname plotGuitar
#'
#' @export

setMethod("plotGuitar",
          "SummarizedExomePeak",
               function(sep,
                        txdb = NULL,
                        save_pdf_prefix = NULL,
                        include_control_regions = TRUE,
                        guitar_coordinate = NULL,
                        save_dir = ".") {
#Check the installation state of the guitar plot
if(!require(Guitar)){
  stop("the 'Guitar' package needs to be installed first from bioconducor.")
}

if(sum(grepl("mod",rownames(sep))) < 10 ) {
  stop("guitar plot function cannot be performed for total peaks number < 10.")
}

stopifnot(!(is.null(txdb) & is.null(guitar_coordinate)))

if(is.null(guitar_coordinate)){
  guitar_coordinate = quiet( Guitar::makeGuitarTxdb(txdb) )
}

if(is.null(save_pdf_prefix)) {
  dir_arg <- NA
} else {
  if(!dir.exists(save_dir)) {
    dir.create(save_dir)
  }
  dir_arg <- file.path( save_dir, save_pdf_prefix)
}

  #first check whether the input object contains quantification result.
    #if so, we need only plot the modification peaks and control peaks.
if(nrow(DESeq2Results(sep)) == 0){
row_grl <- rowRanges( sep )

gr_list <- list(
  peaks = row_grl[grepl("mod",names(row_grl) )],
  background = row_grl[grepl("control",names(row_grl) )]
)

if(!include_control_regions){
  gr_list <- gr_list[-2]
}

quiet(
suppressWarnings(
              GuitarPlot(
                 gfeatures = gr_list,
                 GuitarCoordsFromTxDb = guitar_coordinate,
                 saveToPDFprefix = dir_arg
                        )
                )
)


} else {
  #In case of the collumn design contains only modification
  if(!any(sep$design_Treatment)){

    row_grl <- rowRanges( sep )

    indx_sig <- which( DESeq2Results(sep)$padj < .05 & DESeq2Results(sep)$log2FoldChange > 0 )

    gr_lab = "mod padj < .05"

    if( length(indx_sig) < floor( sum(grepl("mod_", rownames(sep))) * 0.01 ) ){

      indx_sig <- which( DESeq2Results(sep)$pvalue < .05 & DESeq2Results(sep)$log2FoldChange > 0 )

      gr_lab = "mod p < .05"

    }

    gr_list <- list(mod_peaks = row_grl[ grepl("mod", rownames(sep)) ][indx_sig],
                    background = row_grl[ grepl("control", names(row_grl)) ]
    )

    names(gr_list)[1] <- gr_lab

    rm(indx_sig, gr_lab)

   if(!include_control_regions){
     gr_list <- gr_list[-2]
   }

  quiet(
   suppressWarnings(

                  Guitar::GuitarPlot(
                  gfeatures = gr_list,
                  GuitarCoordsFromTxDb = guitar_coordinate,
                  saveToPDFprefix = dir_arg
                )

            )
  )


  } else {
    #Finally, the guitar plot will be generated for differential modification results.

    row_grl <- rowRanges( sep )

    indx_hyper <- which( DESeq2Results(sep)$padj < .05 &
                           DESeq2Results(sep)$log2FoldChange > 0)

    indx_hypo <- which( DESeq2Results(sep)$padj < .05 &
                          DESeq2Results(sep)$log2FoldChange < 0)

    list_names <- c("hyper padj < .05", "hypo padj < .05")

    min_positive <- floor(sum(grepl("mod_", rownames(sep))) * 0.1)

    if(length(indx_hyper) + length(indx_hypo) < min_positive){

      indx_hyper <- which( DESeq2Results(sep)$pvalue < .05 &
                             DESeq2Results(sep)$log2FoldChange > 0)

      indx_hypo <- which( DESeq2Results(sep)$pvalue < .05 &
                            DESeq2Results(sep)$log2FoldChange < 0)

      list_names <- c("hyper p < .05", "hypo p < .05")

    }

    gr_list <- list(hyperMod = row_grl[grepl("mod",rownames(sep))][indx_hyper],
                    hypoMod = row_grl[grepl("mod",rownames(sep))][indx_hypo]
    )

    names(gr_list) <- list_names

    suppressWarnings(
      quiet(
      GuitarPlot(
        gfeatures = gr_list,
        GuitarCoordsFromTxDb = guitar_coordinate,
        saveToPDFprefix = dir_arg
      )
      )
    )

}

}

})
