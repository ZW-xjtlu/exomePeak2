#' @title Plot the Distribution of the Length of Exons Overlapped by the RNA Modification Peaks/Sites
#'
#' @description This function plot the distribution of the exon length for peaks containing exons.
#' @details
#' If the SummarizedExomePeaks object contains LFC statistics, the significantly modified peaks
#' with IP to input log2FC > 0 and DESeq2 Wald test padj < .05 will be plotted .
#'
#' If the SummarizedExomePeaks object contains interactive LFC statistics, both the hyper modification
#' and hypo modification peaks with DESeq2 Wald test p values < .05 will be plotted.
#'
#' @param sep a \code{\link{SummarizedExomePeak}} object.
#' @param txdb a \code{\link{TxDb}} object containing the transcript annotation.
#' @param save_pdf_prefix a \code{character} if provided, a pdf file with the given name will be saved under the current working directory.
#' @param include_control_regions a \code{logical} for whether to include the control regions or not; Default \code{= TRUE}.
#' @param save_dir a \code{character} for the directory to save the plot; Default \code{= "."}.
#'
#' @return a ggplot object
#'
#' @import SummarizedExperiment
#'
#' @aliases plotExonLength
#'
#' @rdname plotExonLength-methods
#'
#' @examples
#'
#'\dontrun{
#'sep #a SummarizedExomePeak object
#'txdb #a TxDb object
#'plotExonLength(sep,txdb)
#'}
#'
#' @export
#'
setMethod("plotExonLength",
          "SummarizedExomePeak",
                  function(sep,
                           txdb = NULL,
                           save_pdf_prefix = NULL,
                           include_control_regions = TRUE,
                           save_dir = ".") {

if( sum(grepl("mod", rownames(sep))) < 10 ) {
    stop("exon length plot cannot be performed for total peaks number < 10.")
}

  stopifnot(!is.null(txdb))

  #first check whether the input object contains any quantification result.

  #if so, we need only plot the modification peaks and control peaks.

  if(is.null(DESeq2Results(sep))){

    row_grl <- rowRanges( sep )

    gr_list <- list(
      peaks = row_grl[grepl("mod",names(row_grl) )],
      control = row_grl[grepl("control",names(row_grl) )]
    )

    if(!include_control_regions){
      gr_list <- gr_list[-2]
    }

    suppressWarnings(

      exonPlot(
        gfeatures = gr_list,
        txdb = txdb,
        save_pdf_prefix = save_pdf_prefix,
        save_dir = save_dir
      )

    )

  } else {
    #In case of the collumn design contains only modification
    if(!any(sep$design_Treatment)) {

      row_grl <- rowRanges( sep )

      indx_sig <- which( DESeq2Results(sep)$padj < .05 & DESeq2Results(sep)$log2FoldChange > 0 )

      gr_lab = "mod padj < .05"

      if( length(indx_sig) < floor( sum(grepl("mod_", rownames(sep))) * 0.01 ) ){

      indx_sig <- which( DESeq2Results(sep)$pvalue < .05 & DESeq2Results(sep)$log2FoldChange > 0 )

      gr_lab = "mod p < .05"

      }

      gr_list <- list(mod_peaks = row_grl[grepl("mod",rownames(sep))][indx_sig],
                      control = row_grl[grepl("control",names(row_grl) )]
      )

      names(gr_list)[1] <- gr_lab

      rm(indx_sig, gr_lab)

      if(!include_control_regions){
        gr_list <- gr_list[-2]
      }

      suppressWarnings(

        exonPlot(
          gfeatures = gr_list,
          txdb = txdb,
          save_pdf_prefix = save_pdf_prefix,
          save_dir = save_dir
        )

      )


    } else {

      row_grl <- rowRanges( sep )

      indx_hyper <- which( DESeq2Results(sep)$padj < .05 &
                             DESeq2Results(sep)$log2FoldChange > 0)

      indx_hypo <- which( DESeq2Results(sep)$padj < .05 &
                            DESeq2Results(sep)$log2FoldChange < 0)

      list_names <- c("hyper padj < 0.05", "hypo padj < 0.05")

      min_positive <- floor(sum(grepl("mod_", rownames(sep))) * 0.1)

      if(length(indx_hyper) + length(indx_hypo) < min_positive){

      indx_hyper <- which( DESeq2Results(sep)$padj < .05 &
                           DESeq2Results(sep)$log2FoldChange > 0)

      indx_hypo <- which( DESeq2Results(sep)$padj < .05 &
                          DESeq2Results(sep)$log2FoldChange < 0)

      list_names <- c("hyper p < 0.05", "hypo p < 0.05")

      }

      gr_list <- list(hyperMod = row_grl[grepl("mod",rownames(sep))][indx_hyper],
                      hypoMod = row_grl[grepl("mod",rownames(sep))][indx_hypo]
      )

      if( any( elementNROWS(gr_list) == 0 ) ){

      gr_list <- c(gr_list,
                   background = row_grl[grepl("control",rownames(sep))])

      }

      names(gr_list)[1:2] <- list_names

      suppressWarnings(

        exonPlot(
          gfeatures = gr_list,
          txdb = txdb,
          save_pdf_prefix = save_pdf_prefix,
          save_dir = save_dir
        )

     )

  }

}
})
