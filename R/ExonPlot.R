#' @title plot the distribution for the length of the overlapped exons.
#' @param gfeatures a list of GRanges or GRangesList.
#' @param txdb a txdb object.
#' @param save_pdf_prefix provided to save a pdf file.
#' @param save_dir a character indicating the directory to save the plot; default ".".
#' @importFrom GenomicFeatures exons
#' @importFrom IRanges subsetByOverlaps
#' @import ggplot2
#'
#'
exonPlot <- function(
  gfeatures,
  txdb,
  save_pdf_prefix = NULL,
  save_dir = "."
) {

  ex_txdb <- exons(txdb)

  list_exlengths <- lapply( gfeatures, function(x) width( subsetByOverlaps( ex_txdb, x) ) )

  plot_df <- data.frame( log2_exon_lengths = log2( unlist(list_exlengths,use.names = FALSE) ),
                         region = rep(names(list_exlengths),sapply(list_exlengths,length) ) )

  plot_df <- plot_df[plot_df$log2_exon_lengths >= 2 &  plot_df$log2_exon_lengths <= 15, ]

  p1 <- ggplot(plot_df,
               aes(x = log2_exon_lengths,
                   fill = region)) +
    geom_density(alpha = .5,
                 linetype = 0) +
    theme_classic() +
    labs(x = "log2 exon lengths", title = "Exon Length Distribution") +
    xlim(2,15) +
    theme(axis.text.x = element_text(face = "bold",
                                     colour = "darkblue"),
          axis.text.y = element_text(face = "bold",
                                     colour = "darkblue"),
          plot.margin = margin(t = 1.5,
                               r = 0.5,
                               b = 1.5,
                               l = 0.5,
                               unit = "cm"),
          plot.title = element_text(face = "bold",
                                    vjust = 1)) +
    scale_fill_brewer(palette = "Dark2")

 if(!is.null(save_pdf_prefix)) {

    if(!dir.exists(save_dir)) {
      dir.create(save_dir)
    }

    ggsave( file.path( save_dir, paste0(save_pdf_prefix, "_exl.pdf")),
                                        p1, width = 5, height = 3)

  }

  return(p1)
}

