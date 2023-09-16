## Plot Poisson GLM fits for GC content bias detection
##
## return save PDF files of GC bias curves
##
plotGCbias <- function(se,
                       fig_dir = NULL,
                       savePrefix = "gc_fit"){
  #require(ggplot2)
  if (!dir.exists(fig_dir)) dir.create(fig_dir)
  plot_df <- data.frame(glmFit = as.numeric(metadata(se)[["fitm"]]),
                        sample = rep(paste0("sample_", seq_len(ncol(se))), each = 200),
                        IP_input = rep(se$IP_input, each = 200),
                        gc = rep(seq(quantile(rowData(se)$gc, 0.05, na.rm = TRUE),
                                     quantile(rowData(se)$gc, 0.95, na.rm = TRUE),
                                     length.out = 200), ncol(se)))

  if(!is.null(se$Perturbation)) {
    plot_df$Perturbation <- rep(se$Perturbation, each = 200)
    ggplot(plot_df, aes(x=gc, y=glmFit, group = sample, colour = IP_input, linetype = Perturbation)) +
      geom_line() + scale_color_manual(values = c('#44AA99', '#332288')) +
      theme_bw() + labs(x = "GC content", y="Normalized Coverage Fit")
  }else{
    ggplot(plot_df, aes(x=gc, y=glmFit, group = sample, colour = IP_input)) +
      geom_line() + scale_color_manual(values = c('#44AA99', '#332288')) +
      theme_bw() + labs(x = "GC content", y="Normalized Coverage Fit")
  }

  ggsave(file.path(fig_dir,paste0(savePrefix,".pdf")), width = 5, height = 3)
}

