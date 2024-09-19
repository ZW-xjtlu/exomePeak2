*exomePeak2* 
================
2021-02-07

**exomePeak2** provides peak detection and differential methylation for Methylated RNA Immunoprecipitation Sequencing (**MeRIP-Seq**) data. *MeRIP-Seq* is a commonly applied sequencing assay that measures the location and abundance of RNA modification sites under specific cellular conditions. In practice, the technique is sensitive to PCR amplification biases commonly found in NGS data. In addition, the efficiency of immunoprecipitation often varies between different IP samples. *exomePeak2* can perform peak calling and differential analysis independent of GC content bias and IP efficiency changes. 

To install exomePeak2 from Github, type the following command in R console.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Rsamtools", "GenomicAlignments", "GenomicRanges", 
                       "GenomicFeatures", "DESeq2", "ggplot2", "mclust", "BSgenome", 
                       "Biostrings", "GenomeInfoDb", "BiocParallel", "IRanges", 
                       "S4Vectors", "rtracklayer", "methods", "stats", 
                       "utils", "BiocGenerics", "magrittr", "speedglm", "splines", "txdbmaker“))

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("ZW-xjtlu/exomePeak2", build_vignettes = TRUE)
```
To view the documentation of exomePeak2, type `browseVignettes("exomePeak2")` after installation.
