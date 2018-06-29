#'@export
plot_size_factors <- function(sep){
 require(ggplot2)
 require(reshape2)
 plot_df <- sapply(c("Control","Both","Methylation"),function(x)estimate_size_factors(sep,from = x)$SE$sizeFactor)
 plot_df <- reshape2::melt(plot_df)
 colnames(plot_df) <- c("bam_files","Estimation_Methods","size_factors")
 plot_df$Estimation_Methods <- factor(plot_df$Estimation_Methods, levels = c("Control","Both","Methylation"))
 plot_df$bam_files <- as.factor( plot_df$bam_files )
 plot_df$IP_input <- "input"
 plot_df$IP_input[rep(sep$SE$design_IP ,3)] <- "IP"
 plot_df$Treatment <- "untreated"
 plot_df$Treatment[rep(sep$SE$design_Treatment ,3)] <- "treated"
 plot_df$samples <- paste0(plot_df$Treatment,"_",plot_df$IP_input)

 Rep_marks <- paste0("Rep", rep(sequence((table(plot_df$samples)/3)[unique(plot_df$samples)]),3))

 plot_df$samples <- paste0(plot_df$samples,"_", Rep_marks)
 ggplot(plot_df,aes(x = samples, y = size_factors)) +
   geom_bar(stat = "identity",
            position = "dodge",
            aes(fill = Estimation_Methods),
            width = 0.8,
            colour = "black") +
   theme_classic() +
  scale_fill_brewer(palette = "Dark2") +
   theme(axis.text.x = element_text(angle = 310,hjust = 0,face = "bold",colour = "darkblue")) +
   labs(x = "Samples", y = "Size factors", title = "Compare size factors estimated by different methods") +
   theme( plot.margin = margin(t = 1,
                               r = 0.5,
                               b = 0.5,
                               l = 1.5,
                               unit = "cm") )

#model_df <- as.data.frame( sapply(c("Control","Both","Methylation"),function(x)estimate_size_factors(sep,from = x)$sizeFactors))
#model_df$Sum_reads <- colSums(assay(sep$SE))
#round( cor(model_df) , 3)
}
