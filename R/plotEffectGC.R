#'@title Plot the relationship between the effect sizes (GLM log2FC estimates) and GC content.
#'@description This function plot the scatter plot between GC content and the methylation / differential methylation levels of the modification sites.
#'
#'@details By default, this function will generate the countour of the scatter plot, and a linear regression line indicating the trend between
#'GC content and log2 fold change or log2 odds ratio returned by DESeq2.
#'The significant changed methylation sites will be grouped and lebeled in different colours.
#'
#'@param bsgenome a \code{\link{BSgenome}} object for the genome sequence, it could be the name of the reference genome recognized by \code{\link{getBSgenom}}.
#'
#'@param txdb a \code{\link{TxDb}} object for the transcript annotation, it could be the name of the reference genome recognized by \code{\link{makeTxDbFromUCSC}}
#'
#'@param save_pdf_prefix a character, if provided, a pdf file with the given name will be saved under the current directory.
#'
#'@param fragment_length the expected fragment length of the sequencing library; Default 100.
#'
#'@param binding_length the expected antibody binding length of IP; Default 25.
#'
#'@param effective_gc whether to calculate the weighted GC content by the probability of reads alignment; default FALSE.
#'
#'@param drop_overlapped_genes whether to mask the overlapped genes in gene annotation, this is meaningful because the GC content estimation is conducted only on exons; Default TRUE.
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
setMethod("plotEffectGC",
          "SummarizedExomePeak",
                function(sep,
                         bsgenome = NULL,
                         txdb = NULL,
                         save_pdf_prefix = NULL,
                         fragment_length = 100,
                         binding_length = 25,
                         effective_GC = FALSE,
                         drop_overlapped_genes = TRUE) {

if(is.null(colData( sep )$sizeFactor)){
    sep <- estimateSeqDepth(sep)
}

if(is.null(DESeq2Results(sep))) {
  if(any(sep$design_Treatment)) {
    sep <- glmDM(sep)
   } else {
    sep <- glmMeth(sep)
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
                          drop_overlapped_genes = drop_overlapped_genes,
                          effective_GC = effective_GC
)

}


Decision <- rep("Insignificant",nrow(DESeq2Results(sep)))

if(!any(sep$design_Treatment)) {

indx_sig <- which( DESeq2Results(sep)$padj < .05 & DESeq2Results(sep)$log2FoldChange > 0 )

if( length(indx_sig) < floor( sum(grepl("meth_", rownames(sep))) * 0.01 ) ){

indx_sig <- which( DESeq2Results(sep)$pvalue < .05 & DESeq2Results(sep)$log2FoldChange > 0 )

Decision[indx_sig] <- "p value < .05"

} else {

Decision[indx_sig] <- "padj < .05"

}

} else {

  if(length(which(DESeq2Results(sep)$padj < .05)) <
     floor(sum(grepl("meth_", rownames(sep))) * 0.005)) {
    Decision[DESeq2Results(sep)$pvalue < .05] <- "p value < .05"

  } else {
    Decision[DESeq2Results(sep)$padj < .05] <- "p adj < .05"

  }

}

GC_content_meth <- elementMetadata(sep)$GC_content[grepl("meth_",rownames(sep))]

na_idx <- is.na( DESeq2Results(sep)$log2FoldChange ) | is.na(GC_content_meth)

plot_df = data.frame(
  Log2FC = DESeq2Results(sep)$log2FoldChange[!na_idx],
  GC_idx = GC_content_meth[!na_idx],
  Label = Decision[!na_idx]
)

if(!any(sep$design_Treatment)) {
  ylabel <- "log2 IP/input Fold Change"
  mtitle <- "GC Content Against log2 Fold Change Estimate"
} else {
  ylabel <- "log2 Risk Ratio"
  mtitle <- "GC Content Against log2 Risk Ratio Estimate"
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
