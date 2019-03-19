library(exomePeak2)
library(Rsamtools)
library(GenomicAlignments)

example_bam <- readGAlignments("/Users/zhenwei/Documents/GitHub/exomePeak2/bam/SRR1182605_sorted.bam")
example_bam <- example_bam[seqnames(example_bam) == "chr21"]
Rsamtools::as
