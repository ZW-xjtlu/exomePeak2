library(testthat)
library(exomePeak2)

test_that( "Check Bams", {

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
                              paired_end = T)

  MeRIP_Seq_Alignment2 <- scanMeripBAM(bam_files =
                                     c("./bam/SRR1182619.bam",
                                       "./bam/SRR1182621.bam",
                                       "./bam/SRR1182623.bam",
                                       "./bam/SRR1182620.bam",
                                       "./bam/SRR1182622.bam",
                                       "./bam/SRR1182624.bam",
                                       "./bam/SRR1182603.bam",
                                       "./bam/SRR1182605.bam",
                                       "./bam/SRR1182604.bam",
                                       "./bam/SRR1182606.bam"),
                                  design_ip = metadata(MeRIP_Seq_Alignment)$design_IP,
                                  design_treatment = metadata(MeRIP_Seq_Alignment)$design_Treatment,
                                  paired_end = T)

  expect_is(MeRIP_Seq_Alignment, "MeripBamFileList" )
  expect_is(Parameter( MeRIP_Seq_Alignment ),"ScanBamParam")
  expect_is(RandomPrimer( MeRIP_Seq_Alignment ),"logical")
  expect_equal(MeRIP_Seq_Alignment, MeRIP_Seq_Alignment2)

} )


