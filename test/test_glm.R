library(testthat)

test_that( "Quantification", {
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(exomePeak2)

  MeRIP_Seq_Alignment <- scan_merip_bams(
    bam_ip = c("./bam/SRR1182619.bam",
               "./bam/SRR1182621.bam",
               "./bam/SRR1182623.bam"),
    bam_input = c("./bam/SRR1182620.bam",
                  "./bam/SRR1182622.bam",
                  "./bam/SRR1182624.bam"),
    bam_treated_ip = c("./bam/SRR1182603.bam",
                       "./bam/SRR1182605.bam"),
    bam_treated_input = c("./bam/SRR1182604.bam",
                          "./bam/SRR1182606.bam"),
    paired_end = TRUE
  )

  SummarizedExomePeaks <- merip_peak_calling(MeRIP_Seq_Alignment,txdb = TxDb.Hsapiens.UCSC.hg19.knownGene)

  SEP_meth <- readRDS("SEP_meth.rds")
  SEP_meth <- estimateSeqDepth(SEP_meth)
  plotSizeFactors(SEP_meth)
  SEP_meth <- GCnormalization(SEP_meth,bsgenome = Hsapiens)
  SEP_meth <- glmMeth(SEP_meth)
  SEP_meth <- glmDM(SEP_meth)

  saveRDS(SEP_meth,"SEP_meth.rds")

 } )
