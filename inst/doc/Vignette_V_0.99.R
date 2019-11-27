## ----para, echo = FALSE, results='hide'------------------------------------
BiocStyle::markdown()
knitr::opts_chunk$set(dev="png",fig.show="hold",
               fig.width=8,fig.height=4.5,fig.align="center",
               message=FALSE,collapse=TRUE)
set.seed(1)

## ---- eval = TRUE----------------------------------------------------------
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

## ---- eval = TRUE----------------------------------------------------------
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

## ---- eval = TRUE----------------------------------------------------------
f2 = system.file("extdata", "mod_annot.rds", package="exomePeak2")

MOD_ANNO_GRANGE <- readRDS(f2)

exomePeak2(bam_ip = IP_BAM,
           bam_input = INPUT_BAM,
           gff_dir = GENE_ANNO_GTF,
           genome = "hg19",
           paired_end = FALSE,
           mod_annot = MOD_ANNO_GRANGE)

## ---- eval = TRUE----------------------------------------------------------
MeRIP_Seq_Alignment <- scanMeripBAM(
                         bam_ip = IP_BAM,
                         bam_input = INPUT_BAM,
                         paired_end = FALSE
                        )

## ---- eval = TRUE----------------------------------------------------------
MeRIP_Seq_Alignment <- scanMeripBAM(
    bam_ip = IP_BAM,
    bam_input = INPUT_BAM,
    bam_treated_input = TREATED_INPUT_BAM,
    bam_treated_ip = TREATED_IP_BAM,
    paired_end = FALSE
  ) 

## ---- eval = TRUE----------------------------------------------------------
SummarizedExomePeaks <- exomePeakCalling(merip_bams = MeRIP_Seq_Alignment,
                                         gff_dir = GENE_ANNO_GTF,
                                         genome = "hg19") 

## ---- eval = TRUE----------------------------------------------------------
SummarizedExomePeaks <- exomePeakCalling(merip_bams = MeRIP_Seq_Alignment,
                                         gff_dir = GENE_ANNO_GTF,
                                         genome = "hg19",
                                         mod_annot = MOD_ANNO_GRANGE) 

## ---- eval = TRUE----------------------------------------------------------
SummarizedExomePeaks <- normalizeGC(SummarizedExomePeaks)

## ---- eval = FALSE---------------------------------------------------------
#  SummarizedExomePeaks <- glmM(SummarizedExomePeaks)

## ---- eval = TRUE----------------------------------------------------------
SummarizedExomePeaks <- glmDM(SummarizedExomePeaks)

## ---- eval = TRUE, fig.align='center', fig.height = 2.8, fig.width = 5-----
plotLfcGC(SummarizedExomePeaks) 

## ---- eval = TRUE----------------------------------------------------------
plotSizeFactors(SummarizedExomePeaks)

## ---- eval = TRUE----------------------------------------------------------
exportResults(SummarizedExomePeaks, format = "BED") 

## --------------------------------------------------------------------------
sessionInfo()

