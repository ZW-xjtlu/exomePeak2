#'@title Plot the relationship between the glm estimate and GC content.
#'@description This function plot the scatter plot between GC content and the methylation / differential methylation levels of the modification sites.
#'
#'@details By default, this function will generate the countour of the scatter plot, and a linear regression line indicating the trend between
#'GC content and log2 fold change or log2 odds ratio returned by DESeq2.
#'The significant changed methylation sites will be grouped and lebeled in different colours.
#'
#'@param sep a \code{summarizedExomePeak} object.
#'@param bsgenome a \code{\link{BSgenome}} object for the genome sequence.
#'@param save_pdf_prefix a character, if provided, a pdf file with the given name will be saved under the current directory.
#'@param fragment_length the expected fragment length of the sequencing library.
#' The widths of the row features for quantifying GC content will be re-sized into the fragment length if it falls bellow it.
#'
#'@return a ggplot object.
#'
#'@import ggplot2
#'@import BSgenome
#'
#'@docType methods
#'
#'@name plotBetaGC
#'
#'@rdname plotBetaGC
#'
#'@export
setMethod("plotBetaGC",
          "SummarizedExomePeak",
                function(sep,
                         bsgenome,
                         save_pdf_prefix = NULL,
                         fragment_length = 100) {

if(is.null(colData( sep )$sizeFactor)){
    sep <- estimateSeqDepth(sep)
}

if(is.null(DESeq2Results(sep))){
  if(any(sep$design_Treatment)){
    sep <- glmDM(sep)
  }else{
    sep <- glmMeth(sep)
  }
}

GC_contents <- GC_content_over_grl(
    grl = rowRanges(sep)[grepl("meth",rownames(sep))],
    bsgenome = bsgenome,
    fragment_length = fragment_length
)

Decision <- rep("Insignificant",nrow(DESeq2Results(sep)))

if(!any(sep$design_Treatment)) {
Decision[DESeq2Results(sep)$padj < 0.05] <- "padj < 0.05"
} else {
Decision[DESeq2Results(sep)$pvalue < 0.05] <- "p < 0.05"
}

na_idx <- is.na( DESeq2Results(sep)$log2FoldChange ) | is.na(GC_contents$GC_content)

plot_df = data.frame(
  Log2FC = DESeq2Results(sep)$log2FoldChange[!na_idx],
  GC_idx = GC_contents$GC_content[!na_idx],
  Label = Decision[!na_idx]
)

if(!any(sep$design_Treatment)) {
  ylabel <- "log2 beta values"
  mtitle <- "GC content against beta estimates"
} else {
  ylabel <- "log2 odds ratios"
  mtitle <- "GC content against log2 odds ratios"
}

plot_df$GC_idx <- as.numeric(plot_df$GC_idx)

plot_df <- plot_df[plot_df$GC_idx < 0.88 & plot_df$GC_idx > 0.2,]


p1 <- ggplot(plot_df, aes(x =  GC_idx , y = Log2FC )) +
                 geom_point(aes(group = Label,
                                colour = Label),
                            size = .05,
                            alpha = .5) +
            theme_classic() +
            scale_colour_manual(values = c("blue", "red")) +
            labs(x = "GC contents",
                 y = ylabel,
                 title = mtitle,
                 subtitle = save_pdf_prefix) +
            xlim(c(0.2,0.9))

#p2 <- ggplot(plot_df, aes(x = GC_idx, fill = Label)) + geom_density(linetype = 0, alpha = .4) + theme_classic() + xlim(c(0.25,0.75)) + scale_fill_brewer(palette = "Dark2") + labs(x = "GC contents", title = "GC content distribution", subtitle = save_pdf_prefix)

#suppressWarnings( suppressMessages( ggsave(paste0(HDER,"_GC_bias.pdf"),p1,width = 5,height = 2.6) ) )
#suppressWarnings( ggsave(paste0(HDER,"_GC_dist.pdf"),p2,width = 5,height = 2.6) )

if(!is.null( save_pdf_prefix )){

suppressMessages( ggsave(
                  paste0(save_pdf_prefix,
                        "_bt_GC.pdf"),
                         p1,
                         width = 4.5,
                         height = 3) )

}

return(p1)

})
