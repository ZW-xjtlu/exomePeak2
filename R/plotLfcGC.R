#'@title Scatter plot between the GLM log2FC estimates and GC content.
#'@description \code{plotLfcGC} plot the scatter plot between GC content and the (differential) modification LFCs.
#'
#'@details By default, this function will generate a scatter plot between GC content and the log2FC value.
#'The significant modification sites will be lebeled in different colours.
#'
#'@param sep a \code{\link{summarizedExomePeak}} object.
#'
#'@param bsgenome a \code{\link{BSgenome}} object for the genome sequence, it could be the name of the reference genome recognized by \code{\link{getBSgenom}}.
#'
#'@param txdb a \code{\link{TxDb}} object for the transcript annotation, it could be the name of the reference genome recognized by \code{\link{makeTxDbFromUCSC}}.
#'
#'@param save_pdf_prefix a \code{character}, if provided, a pdf file with the given name will be saved under the current directory; Default \code{= NULL}.
#'
#'@param fragment_length a \code{numeric} value for the expected fragment length in the RNA-seq library; Default \code{= 100}.
#'
#'@param binding_length a \code{numeric} value for the expected antibody binding length in IP samples; Default \code{= 25}.
#'
#'@param effective_gc a \code{logical} value of whether to calculate the weighted GC content by the probability of reads alignment; default \code{= FALSE}.
#'
#'@param save_dir a \code{character} for the directory to save the plot; default ".".
#'
#'@return a \code{ggplot} object.
#'
#'@examples
#'
#'\dontrun{
#'sep #a SummarizedExomePeak object
#'plotLfcGC(sep)
#'}
#'
#'@import ggplot2
#'@import BSgenome
#'
#'@docType methods
#'
#'@name plotLfcGC
#'@rdname plotLfcGC
#'
#'@export
setMethod("plotLfcGC",
          "SummarizedExomePeak",
                function(sep,
                         bsgenome = NULL,
                         txdb = NULL,
                         save_pdf_prefix = NULL,
                         fragment_length = 100,
                         binding_length = 25,
                         effective_GC = FALSE,
                         save_dir = ".") {

if(is.null(colData( sep )$sizeFactor)){
    sep <- estimateSeqDepth(sep)
}

if(is.null(DESeq2Results(sep))) {
  if(any(sep$design_Treatment)) {
    sep <- glmDM(sep)
   } else {
    sep <- glmM(sep)
  }
}

if(any(is.null(elementMetadata( sep )$GC_content),
       is.null(elementMetadata( sep )$feature_length))) {

stopifnot(!is.null(bsgenome))

stopifnot(!is.null(txdb))

bsgenome <- getBSgenome( bsgenome )

elementMetadata( sep ) <- GC_content_over_grl(
                          bsgenome = bsgenome,
                          txdb = txdb,
                          grl = rowRanges( sep ),
                          fragment_length = fragment_length,
                          binding_length = binding_length,
                          effective_GC = effective_GC
)

}


Decision <- rep("Insignificant",nrow(DESeq2Results(sep)))

if(!any(sep$design_Treatment)) {

indx_sig <- which( DESeq2Results(sep)$padj < .05 & DESeq2Results(sep)$log2FoldChange > 0 )

if( length(indx_sig) < floor( sum(grepl("mod_", rownames(sep))) * 0.01 ) ){

indx_sig <- which( DESeq2Results(sep)$pvalue < .05 & DESeq2Results(sep)$log2FoldChange > 0 )

Decision[indx_sig] <- "p value < 0.05"

} else {

Decision[indx_sig] <- "padj < 0.05"

}

} else {

  if(length(which(DESeq2Results(sep)$padj < .05)) <
     floor(sum(grepl("mod_", rownames(sep))) * 0.1)) {
    Decision[DESeq2Results(sep)$pvalue < .05] <- "p value < 0.05"

  } else {
    Decision[DESeq2Results(sep)$padj < .05] <- "p adj < 0.05"

  }

}

GC_content_mod <- elementMetadata(sep)$GC_content[grepl("mod_",rownames(sep))]

na_idx <- is.na( DESeq2Results(sep)$log2FoldChange ) | is.na(GC_content_mod)

plot_df = data.frame(
  Log2FC = DESeq2Results(sep)$log2FoldChange[!na_idx],
  GC_idx = GC_content_mod[!na_idx],
  Label = Decision[!na_idx]
)

if(!any(sep$design_Treatment)) {
  ylabel <- "IP/input log2 Fold Change"
  mtitle <- "GC Content Against log2 Fold Change Estimates"
} else {
  ylabel <- "Differential log2 Fold Cange"
  mtitle <- "GC Content Against log2 Fold Change Estimates"
}

plot_df$GC_idx <- as.numeric(plot_df$GC_idx)

plot_df <- plot_df[plot_df$GC_idx < 0.88 & plot_df$GC_idx > 0.2,]

p1 <- ggplot(plot_df, aes(x = GC_idx , y = Log2FC )) +
                 geom_point(aes(group = Label,
                                colour = Label),
                            size = .05,
                            alpha = .5) +
            theme_classic() +
            scale_colour_manual(values = c("blue", "red")) +
            labs(x = "GC content",
                 y = ylabel,
                 title = mtitle,
                 subtitle = save_pdf_prefix) +
            xlim(c(0.2,0.9))

if(!is.null( save_pdf_prefix )) {

  if(!dir.exists(save_dir)) {
    dir.create(save_dir)
  }

ggsave(file.path(save_dir,
                 paste0(save_pdf_prefix, "_lfc_GC.pdf")),
                 p1,
                 width = 4.5,
                height = 3)

}

return(p1)

})
