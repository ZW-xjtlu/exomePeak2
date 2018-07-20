counts_peaks <- assay(sep)[grepl("meth",rownames(sep)),]
counts_peaks <- counts_peaks[rowSums(counts_peaks) > 10,]
system.time(
Test <- apply(counts_peaks,
      1,
      function(x){
        summary(glm(x~sep$design_IP, family = poisson(link = "log")))$coefficients[2,4]
      })
)


raw_bins <- readRDS("Raw_Count_bins.rds")

raw_bins <- raw_bins[rowSums(assay(raw_bins)) > 10,]

colData(raw_bins) <- colData(sep)

Cov = ~ design_IP

dds = suppressMessages( DESeqDataSet(se = raw_bins, design = Cov) )

system.time(
dds <- suppressMessages( DESeq(dds) )
)

#355.336 seconds for the test on 0.8 million bins
dds_results <- results(dds)
#kind of linear, and DESeq2 only cost maximum 10 min.

#So, peak calling should be conducted with GC correction + DESeq2.

