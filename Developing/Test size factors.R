#Investigate of the optimal size factor estimation method in MeRIP-seq data
library(SummarizedExperiment)
library(ggplot2)

SE_Peak_counts <- readRDS("./Developing/Bins_100_10.rds")
colData(SE_Peak_counts) <- colData( readRDS("./Developing/SEP_dm.rds")$SE )

#The reason it is hard for the CDF method is that the size factors need to be calculated sample wised.
#The calculation method using CDF requires the IP and input being paired, and the pairing is very hard to be defined.
#What if the sum of all IP v.s. sum of all input? brilliant.

IPenrich_graph <- function(SE_peaks,resl_plot = 1000) {
stopifnot(any(SE_peaks$design_IP)&any(!SE_peaks$design_IP))

rowgr <- unlist(rowRanges(SE_peaks))

gene_id_indx <- rowgr$gene_id[cumsum( elementNROWS(rowRanges(SE_peaks)) )]

rm(rowgr)

if(sum(SE_peaks$design_IP) > 1) {
  IP_total <- rowSums( assay(SE_peaks)[,SE_peaks$design_IP] )
} else {
  IP_total <- as.vector( assay(SE_peaks)[,SE_peaks$design_IP] )
}

if(sum(!SE_peaks$design_IP) > 1) {
input_total <- rowSums(assay(SE_peaks)[,!SE_peaks$design_IP])
} else {
input_total <- as.vector( assay(SE_peaks)[,!SE_peaks$design_IP] )
}

exp_IP <- tapply(IP_total,gene_id_indx,sum)
exp_input <- tapply(input_total,gene_id_indx,sum)

IP_total <- IP_total/(rep(exp_IP,table(gene_id_indx)) + 1)
input_total <- input_total/(rep(exp_input,table(gene_id_indx)) + 1)

signal_bg_cutoff <- which.max( cumsum( sort(input_total) ) - cumsum( sort(IP_total) ) )

IP_signal_index <- IP_total > IP_total[ names(sort(IP_total))[signal_bg_cutoff] ]
input_signal_index <- input_total > input_total[ names(sort(input_total))[signal_bg_cutoff] ]

#Sparsing the data
sparse_index <- seq(1,nrow(SE_peaks),length.out = resl_plot)

IP_total <- cumsum(sort(IP_total))[sparse_index]
input_total <- cumsum(sort(input_total))[sparse_index]

Plot_df <- data.frame(
    y = c(IP_total/max(IP_total),input_total/max(input_total)),
    group = rep(c("IP","input"),each = resl_plot),
    x = rep( (1:resl_plot)/resl_plot, 2 )
)

Plot_df$cut_off = Plot_df$x[ which.max( (input_total/max(input_total)) - (IP_total/max(IP_total)) ) ]

p1 <- ggplot(Plot_df) +
  geom_path(aes(x = x, y = y, colour = group), size = 0.8) +
  geom_vline(aes(xintercept = cut_off), size = 0.8, linetype = 2, colour = "green") +
  theme_classic() +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Percentage of bins",
      y = "Percentage of normalized tags",
      title = "Culmulative percentage enrichment in each channel")

ggsave("cumulative_enrichment_plot.pdf",
       p1,
       width = 5,
       height = 3)

return(data.frame(
  IP_signal_index = IP_signal_index,
  input_signal_index = input_signal_index
))

}

indexes <- IPenrich_graph(SE_Peak_counts)

#Generate size factors on those ranges

size_factors_cdf <- estimateSizeFactorsForMatrix( assay(SE_Peak_counts)[!indexes$IP_signal_index,] )


#Check them through guitar

bins_IP_signal <- reduce( unlist( rowRanges(SE_Peak_counts)[indexes$IP_signal_index] ) )
bins_IP_background <- reduce( unlist( rowRanges(SE_Peak_counts)[!indexes$IP_signal_index] ) )

length(bins_IP_signal) #Very reasonable signal regions in number....

Guitar::GuitarPlot( list(
  bins_IP_signal = bins_IP_signal,
  bins_IP_background = bins_IP_background
), readRDS("/Users/zhenwei/Datasets/Gtcoords/Gtcoord_hg19.rds"),saveToPDFprefix = "CDF_method")

#Turns out to be effective?
sum(width( bins_IP_signal))
#cover 60346152 bp.
#cover only 6673872 bp after DESeq call.

sum(width( reduce(unlist(rowRanges(SEP_meth$SE)[grepl("meth", rownames( SEP_meth$SE ) )] ))))

SEP_meth <- readRDS("./Developing/SEP_meth.rds")

#How about culmulative method + DESeq2 + CQN directly?

table( reduce(unlist(rowRanges(SEP_meth$SE)[grepl("meth", rownames( SEP_meth$SE ) )] )) %over% bins_IP_signal )

#Wierd thing: there are 4788 sites not mapped to signals for no reason.....

Missed_sites <- subsetByOverlaps(reduce(unlist(rowRanges(SEP_meth$SE)[grepl("meth", rownames( SEP_meth$SE ) )] )),bins_IP_signal,invert = TRUE)

Guitar::GuitarPlot( list(
  Missed_sites = Missed_sites
), readRDS("/Users/zhenwei/Datasets/Gtcoords/Gtcoord_hg19.rds"),saveToPDFprefix = "Missed_sites")

#It turns out to be the sites predominantly on 5'UTR.

##Test this approach on some other data

SE_Peak_counts <- readRDS("./Developing/Bins_100_10.rds")
colData(SE_Peak_counts) <- colData( readRDS("./Developing/SEP_dm.rds")$SE )

indexes <- IPenrich_graph(SE_Peak_counts)

bins_IP_signal <- reduce( unlist( rowRanges(SE_Peak_counts)[indexes$IP_signal_index] ) )

length(bins_IP_signal) #Moderately more signal regions in number.... 85155

SE_Peak_counts  <- SE_Peak_counts[,c(1,4)]

indexes <- IPenrich_graph(SE_Peak_counts)

bins_IP_signal <- reduce( unlist( rowRanges(SE_Peak_counts)[indexes$IP_signal_index] ) )

length(bins_IP_signal) #significant more signal regions in number.... 161944

#How about set a cutoff, if the index bigger than 200000, then it will re-subset to the top 200000 bins.

#Top N local/whole gene enrichment scores.




