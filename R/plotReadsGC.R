#' @title Plot the relationship between GC content and reads abundance.
#'
#' @description  This function plot the local regression curves of the normalized feature abundance against the local GC content levels.
#'
#' @details The read abundances of both the control and the modification site regions are plotted,
#' the read counts are normalized using the following method:
#'
#'  normalized feature abundance = ( ( read_count / size_factor) / region_length ) * 500
#'
#'  By default, it will use the sequencing depth size factor defined in \code{summarizedExomePeak},
#'  if the size factor is not found, new size factors will be estimated using the default method in \code{estimateSeqDepth}.
#'
#'@param bsgenome a \code{\link{BSgenome}} object for the genome sequence, it could be the name of the reference genome recognized by \code{\link{getBSgenom}}.
#'
#'@param txdb a \code{\link{TxDb}} object for the transcript annotation, it could be the name of the reference genome recognized by \code{\link{makeTxDbFromUCSC}}.
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
#'@param pool_replicates a logical indicating whether to pool the replicates in the local regression fit; default FALSE.
#'
#'@return a ggplot object
#'
#'@import ggplot2
#'@import BSgenome
#'@import SummarizedExperiment
#'@importFrom reshape2 melt
#'
#'@docType methods
#'
#'@name plotReadsGC
#'
#'@rdname plotReadsGC
#'
#'@export
setMethod("plotReadsGC",
          "SummarizedExomePeak",
               function(sep,
                        bsgenome = NULL,
                        txdb = NULL,
                        save_pdf_prefix = NULL,
                        fragment_length = 100,
                        binding_length = 25,
                        effective_GC = FALSE,
                        drop_overlapped_genes = TRUE,
                        pool_replicates = FALSE) {

if(is.null(colData( sep )$sizeFactor)){
    sep <- estimateSeqDepth( sep )
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

normalized_counts <- ( t( t(assay(sep))/sep$sizeFactor ) / elementMetadata( sep )$feature_length )

if(!is.null(GCsizeFactors( sep ))) {

cqnormalized_counts <- assay(sep)/exp( GCsizeFactors( sep ) )

}

#Remove the rows with average < 50

elementMetadata( sep )$GC_content[rowMeans(assay(sep)) < 50] = NA

Plot_df <- melt(normalized_counts)

Plot_df$GC_cont = rep( elementMetadata( sep )$GC_content, ncol(sep) )

IP_input <- rep( "input", ncol(sep) )

IP_input[ sep$design_IP ] <- "IP"

if(any(sep$design_Treatment)) {

Treatment <- rep( "control", ncol(sep) )

Treatment[ sep$design_Treatment ] <- "treated"

samples <- paste0(Treatment, "_", IP_input)

} else {

samples <- IP_input

}

Rep_marks <- paste0("Rep", sequence((table(samples))[unique(samples)]))

samples <- paste0(samples,"_", Rep_marks)

Plot_df$samples = rep(samples, each = nrow(sep))

Plot_df$IP_input = rep(IP_input, each = nrow(sep))

if(any(sep$design_Treatment)) {
Plot_df$Treatment = rep(Treatment, each = nrow(sep))
}

group <- rep("background",nrow(sep))

group[grepl("meth",rownames( sep ))] <- "methylated"

Plot_df$group <- rep(group,  ncol(sep))

if(!is.null(GCsizeFactors( sep ))){
Plot_df_cqn = melt(cqnormalized_counts)
Plot_df_cqn = cbind(Plot_df_cqn,Plot_df[,4:ncol(Plot_df)])
Plot_df$norm = "original"
Plot_df_cqn$norm = "gc_corrected"
Plot_df = rbind(Plot_df,Plot_df_cqn)
Plot_df$norm = factor(Plot_df$norm,levels = c("original","gc_corrected"))
}

Plot_df_temp <- Plot_df

Plot_df_temp$group <- "all"

Plot_df <- rbind(Plot_df,
                 Plot_df_temp)

rm(Plot_df_temp)

Plot_df <- na.omit(Plot_df)

Plot_df <- Plot_df[Plot_df$GC_cont >= 0.21 & Plot_df$GC_cont <= 0.85,]

if(pool_replicates){

p1 <- ggplot(Plot_df,
             aes(x = GC_cont,
                 y = value,
                 colour = Treatment,
                 linetype = IP_input)) +
  geom_smooth(alpha = .09,size = .9) +
  theme_bw() +
  labs(x = "GC content",
       y = "normalized feature abundance",
       title = "GC content diagnosis plot")

} else {

p1 <- ggplot(Plot_df,
               aes(x = GC_cont,
                   y = value,
                   colour = samples,
                   linetype = IP_input) ) +
  geom_smooth(alpha = .15,size = .9) +
  theme_bw() +
  labs(x = "GC content",
       y = "normalized feature abundance",
       title = "GC content diagnosis plot")

}

if (ncol(sep) <= 11)  p1 = p1 + scale_color_brewer(palette = "Spectral")

if(!is.null(GCsizeFactors( sep ))) {
  p1 = p1 + facet_grid(norm~group,scales = "free_y")
} else {
  p1 = p1 + facet_grid(~group,scales = "free_y")
}

figheight = 6.55 + .2 * round(ncol( sep )/4)

if(!is.null(GCsizeFactors( sep ))) figwidth = 12

if(!is.null(GCsizeFactors( sep )) & !any(sep$design_Treatment)) {
p1 = p1 +
  theme(plot.margin = margin(1,2,1,2,"cm"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "vertical",
          legend.justification = "center")

} else {

p1 = p1 +
  theme(plot.margin = margin(1,5,1,5,"cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.justification = "center")

}

if(!is.null(save_pdf_prefix)) {

suppressMessages( ggsave(paste0(save_pdf_prefix,"_ab_GC.pdf"), p1, width = figwidth, height = figheight) )

}

suppressMessages( return(p1) )

})
