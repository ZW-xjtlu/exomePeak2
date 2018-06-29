#' @title Plot the relationship between GC content and feature abundance.
#'
#' @description  This function plot the local regression curves of the normalized feature abundance (reads count) against the local GC content levels.
#'
#' @details The read abundances of both the control and the modification site regions are plotted,
#' the read counts are normalized using the following method:
#'
#'  normalized feature abundance = ( ( read_count / size_factor) / region_length ) * 500
#'
#'  By default, it will use the sequencing depth size factor defined in \code{summarizedExomePeak},
#'  if the size factor is not found, new size factors will be estimated using the default method in \code{estimate_size_factors}.
#'
#'@param sep a \code{summarizedExomePeak} object.
#'@param bsgenome a \code{\link{BSgenome}} object for the genome sequence.
#'@param fragment_length the expected fragment length of the sequencing library.
#' The widths of the row features for quantifying GC content will be re-sized into the fragment length if it falls bellow it.
#'@param save_pdf_prefix a character, if provided, a pdf file with the given name will be saved under the current directory.
#'@param combine_replicates a logical indicating whether to pool the replicates in the local regression fit; default FALSE.
#'
#'@return a ggplot object
#'
#'@import ggplot2
#'@import BSgenome
#'@import SummarizedExperiment
#'@importFrom reshape2 melt
#'@export
abundance_GC_plot <- function(sep,
                          bsgenome,
                          fragment_length = 100,
                          save_pdf_prefix = NULL,
                          combine_replicates = FALSE) {

if(is.null(colData( sep$SE )$sizeFactor)){
    sep <- estimate_size_factors(sep)
}

GC_contents <- GC_content_over_grl(
                grl = rowRanges(sep$SE),
                bsgenome = bsgenome,
                fragment_length = fragment_length
               )

normalized_counts <- ( t( t(assay(sep$SE))/sep$SE$sizeFactor ) / GC_contents$Indx_length ) * 500

if(!is.null(sep$FS_sizeFactor)) {

cqnormalized_counts <- assay(sep$SE)/sep$FS_sizeFactor

}

Plot_df <- melt(normalized_counts)

Plot_df$GC_cont = rep( GC_contents$GC_content , ncol(sep$SE) )

IP_input <- rep( "input", ncol(sep$SE) )

IP_input[ sep$SE$design_IP ] <- "IP"

if(any(sep$SE$design_Treatment)) {

Treatment <- rep( "control", ncol(sep$SE) )

Treatment[ sep$SE$design_Treatment ] <- "treated"

samples <- paste0(Treatment, "_", IP_input)
} else {
samples <- IP_input
}

Rep_marks <- paste0("Rep", sequence((table(samples))[unique(samples)]))

samples <- paste0(samples,"_", Rep_marks)

Plot_df$samples = rep(samples, each = nrow(sep$SE))

Plot_df$IP_input = rep(IP_input, each = nrow(sep$SE))

if(any(sep$SE$design_Treatment)) {
Plot_df$Treatment = rep(Treatment, each = nrow(sep$SE))
}

group <- rep("background",nrow(sep$SE))

group[grepl("meth",rownames( sep$SE ))] <- "methylated"

Plot_df$group <- rep(group,  ncol(sep$SE))

if(!is.null(sep$FS_sizeFactor)){
Plot_df_cqn = melt(cqnormalized_counts)
Plot_df_cqn = cbind(Plot_df_cqn,Plot_df[,4:ncol(Plot_df)])
Plot_df$norm = "original"
Plot_df_cqn$norm = "cqn"
Plot_df = rbind(Plot_df,Plot_df_cqn)
Plot_df$norm = factor(Plot_df$norm,levels = c("original","cqn"))
}

Plot_df_temp <- Plot_df

Plot_df_temp$group <- "all"

Plot_df <- rbind(Plot_df,
                 Plot_df_temp)

rm(Plot_df_temp)

Plot_df <- na.omit(Plot_df)

Plot_df <- Plot_df[Plot_df$GC_cont >= 0.21 & Plot_df$GC_cont <= 0.85,]

if(combine_replicates){

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

}else{

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

if (ncol(sep$SE) <= 11)  p1 = p1 + scale_color_brewer(palette = "Spectral")

if(!is.null(sep$FS_sizeFactor)) {
  p1 = p1 + facet_grid(group~norm,scales = "free_y")
} else {
  p1 = p1 + facet_grid(~group,scales = "free_y")
}

figheight = 6.55 + .2 * round(ncol( sep$SE )/4)

if(!is.null(sep$FS_sizeFactor)) figwidth = 8.8

else figwidth = 7.5

if(!is.null(sep$FS_sizeFactor) & !any(sep$SE$design_Treatment)) {
p1 = p1 +
  theme(plot.margin = margin(1,2,1,2,"cm"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "vertical",
          legend.justification = "center")

figwidth = 6.3

}else {

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

}
