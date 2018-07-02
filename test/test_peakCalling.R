library(testthat)

test_that( "Peak Calling", {
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(exomePeak2)

  MeRIP_Seq_Alignment <- scanMeripBAM(
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

  SummarizedExomePeaks <- exomePeakCalling(merip_bams = MeRIP_Seq_Alignment,
                                           txdb = TxDb.Hsapiens.UCSC.hg19.knownGene)

  expect_is(SummarizedExomePeaks,"SummarizedExomePeak")
  saveRDS(SummarizedExomePeaks,"SEP_meth.rds")

  #Test the accessors
  #SummarizedExomePeaks <- readRDS("SEP_meth.rds")
  GCsizeFactors(SummarizedExomePeaks) <- matrix(1:10)
  expect_equal(GCsizeFactors(SummarizedExomePeaks),matrix(1:10))


  MeRIP_Seq_Alignment_meth <- scan_merip_bams(
    bam_ip = c("./bam/SRR1182619.bam",
               "./bam/SRR1182621.bam",
               "./bam/SRR1182623.bam"),
    bam_input = c("./bam/SRR1182620.bam",
                  "./bam/SRR1182622.bam",
                  "./bam/SRR1182624.bam"),
    paired_end = TRUE
  )

  SummarizedExomePeaks_meth <- merip_peak_calling(MeRIP_Seq_Alignment_meth,txdb = TxDb.Hsapiens.UCSC.hg19.knownGene)

  expect_that(SummarizedExomePeaks, is_a("list") ) #Later when we use OOP, we should change the list into our own S4 object.
} )
