#' @title Plot the linear relationships between GC content and reads abundance.
#'
#' @description  \code{plotReadsGC} visualizes the local regression curves between the normalized reads abundance and the local GC content.
#'
#' @details The read abundances of both the control and the modification site regions are plotted,
#' the read counts are normalized using the following method:
#'
#' \deqn{normalized feature abundance = ( ( read count / size factor) / region length ) * 500}
#'
#'  By default, it will use the sequencing depth size factor defined in the \code{\link{SummarizedExomePeak}} object,
#'  if the sequencing depth size factor is not found, new size factors will be estimated with the default method in \code{\link{estimateSeqDepth}}.
#'
#'@param sep a \code{\link{SummarizedExomePeak}} object.
#'
#'@param bsgenome a \code{\link{BSgenome}} object for the genome sequence, it could be the name of the reference genome recognized by \code{\link{getBSgenome}}.
#'
#'@param txdb a \code{\link{TxDb}} object for the transcript annotation, it could be the name of the reference genome recognized by \code{\link{makeTxDbFromUCSC}}.
#'
#'@param save_pdf_prefix a \code{character}, if provided, a pdf file with the given name will be saved under the current directory; Default \code{= NULL}.
#'
#'@param fragment_length a \code{numeric} value for the expected fragment length in the RNA-seq library; Default \code{= 100}.
#'
#'@param binding_length a \code{numeric} value for the expected antibody binding length in IP samples; Default \code{= 25}.
#'
#'@param effective_GC a \code{logical} value of whether to calculate the weighted GC content by the probability of reads alignment; default \code{= FALSE}.
#'
#'@param pool_replicates a \code{logical} value of whether to pool the replicates in the local regression fit; default \code{= FALSE}.
#'
#'@param save_dir a \code{character} for the directory to save the plot; default ".".
#'
#'@return a \code{ggplot} object
#'
#'@examples
#'
#'### Load the example SummarizedExomPeak object
#' f1 = system.file("extdata", "sep_ex.rds", package="exomePeak2")
#'
#' sep <- readRDS(f1)
#'
#' ### Visualize the linear relationships between GC content and normalized reads count
#' plotReadsGC(sep)
#'
#'@import ggplot2
#'@import BSgenome
#'@import SummarizedExperiment
#'@import quantreg
#'@importFrom reshape2 melt
#'
#'@aliases plotReadsGC
#'
#'@rdname plotReadsGC-methods
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
                        pool_replicates = FALSE,
                        save_dir = ".") {

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
    effective_GC = effective_GC
  )

}

#Select the rows with average < 50

sep <- sep[rowMeans(assay(sep)) >= 50,]

normalized_counts <- ( t( t(assay(sep))/sep$sizeFactor ) / elementMetadata( sep )$feature_length )

normalized_counts <- scale(normalized_counts)

if(!is.null(GCsizeFactors( sep ))) {

cqnormalized_counts <- assay(sep) / exp( GCsizeFactors( sep ) )

cqnormalized_counts <- scale(  cqnormalized_counts )

}

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

group[grepl("mod",rownames( sep ))] <- "modification"

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
  p1 = p1 + facet_grid(norm~group)
} else {
  p1 = p1 + facet_grid(~group)
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

  if(!dir.exists(save_dir)) {
    dir.create(save_dir)
  }

  suppressMessages( ggsave(file.path(save_dir, paste0(save_pdf_prefix, "ReadsGC.pdf")),
                           p1, width = figwidth, height = figheight) )

}

suppressMessages( return(p1) )

})
