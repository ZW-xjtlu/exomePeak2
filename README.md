The *exomePeak2* user's guide
================
2019-12-20

Introduction
============

**exomePeak2** provides technically independent quantification and peak detection on Methylated RNA immunoprecipitation sequencing data (**MeRIP-Seq**). *MeRIP-Seq* is a commonly applied sequencing technology to measure the transcriptome-wide location and abundance of RNA modification sites under a given cellular condition. However, the quantification and peak calling in *MeRIP-Seq* are sensitive to PCR amplification bias which is prevalent in next generation sequencing (NGS) techniques. In addition, the **RNA-Seq** based count data exhibits biological variation in small reads count. *exomePeak2* collectively address these challanges by introducing a rich set of robust data science models tailored for MeRIP-Seq. With *exomePeak2*, users can perform peak calling, modification site quantification, and differential analysis with a straightforward one-step function. Alternatively, users could define personalized methods for their own analysis through multi-step functions and diagnostic plots.

Package Installation
====================

To install exomePeak2 from Github, use the following codes.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SummarizedExperiment","cqn","Rsamtools",
                       "GenomicAlignments","GenomicRanges","GenomicFeatures",
                       "DESeq2","ggplot2","mclust",
                       "genefilter","BSgenome","BiocParallel",
                       "IRanges","S4Vectors","quantreg",
                       "reshape2","rtracklayer","apeglm","RMariaDB"))

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("ZW-xjtlu/exomePeak2")
```

Peak Calling
============

For peak calling of *MeRIP-Seq* experiment, exomePeak2 demands the reads alignment results in **BAM** files. Users can specify the biological replicates of the IP and input samples by a character vector of the corresponding **BAM** directories at the arguments `bam_ip` and `bam_input` separately.

In the following example, the transcript annotation is provided using GFF files. Transcript annotation can also be provided by the `TxDb` object. exomePeak2 will automatically download the TxDb if the `genome` argument is filled with the corresponding UCSC genome name.

The genome sequence is required to conduct GC content bias correction. If the `genome` argument is missing ( `= NULL` ), exomPeak2 will perform peak calling without correcting the GC content bias.

``` r
library(exomePeak2)

GENE_ANNO_GTF = system.file("extdata", "example.gtf", package="exomePeak2")

f1 = system.file("extdata", "IP1.bam", package="exomePeak2")
f2 = system.file("extdata", "IP2.bam", package="exomePeak2")
f3 = system.file("extdata", "IP3.bam", package="exomePeak2")
f4 = system.file("extdata", "IP4.bam", package="exomePeak2")
IP_BAM = c(f1,f2,f3,f4)

f1 = system.file("extdata", "Input1.bam", package="exomePeak2")
f2 = system.file("extdata", "Input2.bam", package="exomePeak2")
f3 = system.file("extdata", "Input3.bam", package="exomePeak2")
INPUT_BAM = c(f1,f2,f3)

exomePeak2(bam_ip = IP_BAM,
           bam_input = INPUT_BAM,
           gff_dir = GENE_ANNO_GTF,
           genome = "hg19",
           paired_end = FALSE)
## class: SummarizedExomePeak 
## dim: 31 7 
## metadata(0):
## assays(2): counts GCsizeFactors
## rownames(31): peak_11 peak_13 ... control_13 control_14
## rowData names(2): GC_content feature_length
## colnames(7): IP1.bam IP2.bam ... Input2.bam Input3.bam
## colData names(3): design_IP design_Treatment sizeFactor
```

exomePeak2 will export the modification peaks in formats of **BED** file and **CSV** table, the data will be saved automatically under a folder named by `exomePeak2_output`.

The modification peak statistics are derived from the *β*<sub>*i*, 1</sub> terms in the following linear regression design.

*l**o**g*2(*Q*<sub>*i*, *j*</sub>)=*β*<sub>*i*, 0</sub> + *β*<sub>*i*, 1</sub>*I*(*ρ*(*j*)=*I**P*)+*t*<sub>*i*, *j*</sub>

*Q*<sub>*i*, *j*</sub> is the expected value of reads abundence of modification *i* under sample *j*. *β*<sub>*i*, 0</sub> is the intercept coefficient, *β*<sub>*i*, 1</sub> is the coefficient for IP/input log2 fold change, *I*(*ρ*(*j*)=*I**P*) is the regression covariate that is the indicator variable for the sample *j* being IP sample. *t*<sub>*i*, *j*</sub> is the regression offset that account for the sequencing depth variation and the GC content biases.

Under the default settings, the linear models fitted are the regularized **GLM (Generalized Linear Model)** of NB developed by **DESeq2**. If one of the IP and input group has no biological replicates, Poisson GLMs will be fitted to the modification peaks.

Explaination over the columns of the exported table:

-   ***chr***: the chromosomal name of the peak.
-   ***chromStart***: the start of the peak on the chromosome.
-   ***chromEnd***: the end of the peak on the chromosome.
-   ***name***: the unique ID of the modification peak.
-   ***score***: the -log2 p value of the peak.
-   ***strand***: the strand of the peak on genome.
-   ***thickStart***: the start position of the peak.
-   ***thickEnd***: the end position of the peak.
-   ***itemRgb***: the column for the RGB encoded color in BED file visualization.
-   ***blockCount***: the block (exon) number within the peak.
-   ***blockSizes***: the widths of blocks.
-   ***blockStarts***: the start positions of blocks.
-   ***geneID***: the gene ID of the peak.
-   ***ReadsCount.input***: the reads count of the input sample.
-   ***ReadsCount.IP***: the reads count of the IP sample.
-   ***log2FoldChange***: the estimates of IP over input log2 fold enrichment (coefficient estimates of *β*<sub>*i*, 1</sub>).
-   ***pvalue***: the Wald test p value on the modification coefficient.
-   ***padj***: the adjusted Wald test p value using BH approach.

Differential Modification Analysis
==================================

The code below could conduct differential modification analysis (Comparison of Two Conditions) on exon regions defined by the transcript annotation.

In differential modification mode, exomePeak2 will first perform Peak calling on exon regions using both the control and treated samples. Then, it will conduct the differential modification analysis on peaks reported from peak calling using an interactive GLM.

``` r
f1 = system.file("extdata", "treated_IP1.bam", package="exomePeak2")
TREATED_IP_BAM = c(f1)
f1 = system.file("extdata", "treated_Input1.bam", package="exomePeak2")
TREATED_INPUT_BAM = c(f1)

exomePeak2(bam_ip = IP_BAM,
           bam_input = INPUT_BAM,
           bam_treated_input = TREATED_INPUT_BAM,
           bam_treated_ip = TREATED_IP_BAM,
           gff_dir = GENE_ANNO_GTF,
           genome = "hg19",
           paired_end = FALSE)
## class: SummarizedExomePeak 
## dim: 23 9 
## metadata(0):
## assays(2): counts GCsizeFactors
## rownames(23): peak_10 peak_11 ... control_5 control_6
## rowData names(2): GC_content feature_length
## colnames(9): IP1.bam IP2.bam ... treated_IP1.bam
##   treated_Input1.bam
## colData names(3): design_IP design_Treatment sizeFactor
```

In differential modification mode, exomePeak2 will export the differential modification peaks in formats of **BED** file and **CSV** table, the data will also be saved automatically under a folder named by `exomePeak2_output`.

The peak statistics in differential modification setting are derived from the interactive coefficient *β*<sub>*i*, 3</sub> in the following regression design of the **NB GLM**:

*l**o**g*2(*Q*<sub>*i*, *j*</sub>)=*β*<sub>*i*, 0</sub> + *β*<sub>*i*, 1</sub>*I*(*ρ*(*j*)=*I**P*)+*β*<sub>*i*, 2</sub>*I*(*ρ*(*j*)=*T**r**e**a**t**m**e**n**t*)+*β*<sub>*i*, 3</sub>*I*(*ρ*(*j*)=*I**P*&*T**r**e**a**t**m**e**n**t*)+*t*<sub>*i*, *j*</sub>

Explaination for the additional table columns:

-   ***ModLog2FC\_control***: the modification log2 fold enrichment in the control condition.
-   ***ModLog2FC\_treated***: the modification log2 fold enrichment in the treatment condition.
-   ***DiffModLog2FC***: the log2 Fold Change estimates of differential modification (coefficient estimates of *β*<sub>*i*, 3</sub>).
-   ***pvalue***: the Wald test p value on the differential modification coefficient.
-   ***padj***: the adjusted Wald test p value using BH approach.

Quantification and Statistical Analysis with Single Based Modification Annotation
=================================================================================

exomePeak2 supports the modification quantification and differential modification analysis on single based modification annotation. The modification sites with single based resolution can provide a more accurate mapping of modification locations compared with the peaks called directly from the MeRIP-seq datasets.

Some of the datasets in epitranscriptomics have single based resolution, e.x. Data generated by the *m6A-CLIP-Seq* or *m6A-miCLIP-Seq* techniques. Reads count on the single based modification sites could provide a more accurate and consistent quantification on *MeRIP-Seq* experiments with single based annotation.

exomePeak2 will automatically initiate the mode of single based modification quantification by providing a sigle based annotation file under the argument `mod_annot`.

The single based annotation information should be provided to the exomePeak2 function in the format of a `GRanges` object.

``` r
f2 = system.file("extdata", "mod_annot.rds", package="exomePeak2")

MOD_ANNO_GRANGE <- readRDS(f2)

exomePeak2(bam_ip = IP_BAM,
           bam_input = INPUT_BAM,
           gff_dir = GENE_ANNO_GTF,
           genome = "hg19",
           paired_end = FALSE,
           mod_annot = MOD_ANNO_GRANGE)
## class: SummarizedExomePeak 
## dim: 171 7 
## metadata(0):
## assays(2): '' GCsizeFactors
## rownames(171): peak_1 peak_2 ... control_83 control_84
## rowData names(2): GC_content feature_length
## colnames(7): IP1.bam IP2.bam ... Input2.bam Input3.bam
## colData names(3): design_IP design_Treatment sizeFactor
```

In this mode, exomePeak2 will export the analysis result also in formats of **BED** file and **CSV** table, while each row of the table corresponds to the sites of the annotation `GRanges`.

Peak Calling and Visualization in Multiple Steps
================================================

The exomePeak2 package can achieve peak calling and peak statistics calulation with multiple functions.

**1. Check the bam files of MeRIP-seq data before peak calling.**

``` r
MeRIP_Seq_Alignment <- scanMeripBAM(
                         bam_ip = IP_BAM,
                         bam_input = INPUT_BAM,
                         paired_end = FALSE
                        )
```

For MeRIP-seq experiment with interactive design (contain control and treatment groups), use the following code.

``` r
MeRIP_Seq_Alignment <- scanMeripBAM(
    bam_ip = IP_BAM,
    bam_input = INPUT_BAM,
    bam_treated_input = TREATED_INPUT_BAM,
    bam_treated_ip = TREATED_IP_BAM,
    paired_end = FALSE
  ) 
```

**2. Conduct peak calling analysis on exons using the provided bam files.**

``` r
SummarizedExomePeaks <- exomePeakCalling(merip_bams = MeRIP_Seq_Alignment,
                                         gff_dir = GENE_ANNO_GTF,
                                         genome = "hg19") 
```

Alternatively, use the following code to quantify MeRIP-seq data on single based modification annotation.

``` r
SummarizedExomePeaks <- exomePeakCalling(merip_bams = MeRIP_Seq_Alignment,
                                         gff_dir = GENE_ANNO_GTF,
                                         genome = "hg19",
                                         mod_annot = MOD_ANNO_GRANGE) 
```

**3. Estimate size factors that are required for GC content bias correction.**

``` r
SummarizedExomePeaks <- normalizeGC(SummarizedExomePeaks)
```

**4. Report the statistics of modification peaks using Generalized Linear Model (GLM).**

``` r
SummarizedExomePeaks <- glmM(SummarizedExomePeaks) 
```

Alternatively, If the treated IP and input bam files are provided, `glmDM` function could be used to conduct differential modification analysis on modification Peaks with interactive GLM.

``` r
SummarizedExomePeaks <- glmDM(SummarizedExomePeaks)
```

**5. Generate the scatter plot between GC content and log2 Fold Change (LFC).**

``` r
plotLfcGC(SummarizedExomePeaks) 
```

**6. Generate the bar plot for the sequencing depth size factors.**

``` r
plotSizeFactors(SummarizedExomePeaks)
```

<img src="inst/figures/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

**7. Export the modification peaks and the peak statistics with user decided format.**

``` r
exportResults(SummarizedExomePeaks, format = "BED") 
```

Contact
=======

Please contact the maintainer of exomePeak2 if you have encountered any problems:

**ZhenWei**: <zhen.wei@xjtlu.edu.cn>

Session Info
============

``` r
sessionInfo()
## R version 3.5.3 (2019-03-11)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS Sierra 10.12.6
## 
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] zh_CN.UTF-8/zh_CN.UTF-8/zh_CN.UTF-8/C/zh_CN.UTF-8/zh_CN.UTF-8
## 
## attached base packages:
##  [1] splines   parallel  stats4    stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] BSgenome.Hsapiens.UCSC.hg19_1.4.0 BSgenome_1.50.0                  
##  [3] rtracklayer_1.42.2                Biostrings_2.50.2                
##  [5] XVector_0.22.0                    exomePeak2_0.99.1                
##  [7] cqn_1.28.1                        quantreg_5.51                    
##  [9] SparseM_1.77                      preprocessCore_1.44.0            
## [11] nor1mix_1.3-0                     mclust_5.4.5                     
## [13] SummarizedExperiment_1.12.0       DelayedArray_0.8.0               
## [15] BiocParallel_1.16.6               matrixStats_0.54.0               
## [17] Biobase_2.42.0                    GenomicRanges_1.34.0             
## [19] GenomeInfoDb_1.18.2               IRanges_2.16.0                   
## [21] S4Vectors_0.20.1                  BiocGenerics_0.28.0              
## [23] BiocStyle_2.10.0                 
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.4-1         htmlTable_1.13.1        
##  [3] base64enc_0.1-3          rstudioapi_0.10         
##  [5] MatrixModels_0.4-1       bit64_0.9-7             
##  [7] AnnotationDbi_1.44.0     apeglm_1.4.2            
##  [9] geneplotter_1.60.0       knitr_1.23              
## [11] zeallot_0.1.0            Formula_1.2-3           
## [13] Rsamtools_1.34.1         annotate_1.60.1         
## [15] cluster_2.1.0            BiocManager_1.30.4      
## [17] compiler_3.5.3           httr_1.4.0              
## [19] backports_1.1.5          assertthat_0.2.1        
## [21] Matrix_1.2-17            lazyeval_0.2.2          
## [23] acepack_1.4.1            htmltools_0.4.0         
## [25] prettyunits_1.0.2        tools_3.5.3             
## [27] coda_0.19-3              gtable_0.3.0            
## [29] glue_1.3.1               GenomeInfoDbData_1.2.0  
## [31] reshape2_1.4.3           dplyr_0.8.3             
## [33] Rcpp_1.0.2               bbmle_1.0.20            
## [35] vctrs_0.2.0              xfun_0.8                
## [37] stringr_1.4.0            XML_3.98-1.20           
## [39] zlibbioc_1.28.0          MASS_7.3-51.4           
## [41] scales_1.0.0             hms_0.5.0               
## [43] RMariaDB_1.0.6           RColorBrewer_1.1-2      
## [45] yaml_2.2.0               memoise_1.1.0           
## [47] gridExtra_2.3            ggplot2_3.2.1           
## [49] emdbook_1.3.11           biomaRt_2.38.0          
## [51] rpart_4.1-15             latticeExtra_0.6-28     
## [53] stringi_1.4.3            RSQLite_2.1.2           
## [55] genefilter_1.64.0        checkmate_1.9.4         
## [57] GenomicFeatures_1.34.8   rlang_0.4.1             
## [59] pkgconfig_2.0.3          bitops_1.0-6            
## [61] evaluate_0.14            lattice_0.20-38         
## [63] purrr_0.3.2              GenomicAlignments_1.18.1
## [65] htmlwidgets_1.5.1        labeling_0.3            
## [67] bit_1.1-14               tidyselect_0.2.5        
## [69] plyr_1.8.4               magrittr_1.5            
## [71] DESeq2_1.22.2            R6_2.4.0                
## [73] Hmisc_4.2-0              DBI_1.0.0               
## [75] pillar_1.4.2             foreign_0.8-71          
## [77] survival_2.44-1.1        RCurl_1.95-4.12         
## [79] nnet_7.3-12              tibble_2.1.3            
## [81] crayon_1.3.4             rmarkdown_1.14          
## [83] progress_1.2.2           locfit_1.5-9.1          
## [85] grid_3.5.3               data.table_1.12.2       
## [87] blob_1.2.0               digest_0.6.22           
## [89] xtable_1.8-4             numDeriv_2016.8-1.1     
## [91] munsell_0.5.0
```
