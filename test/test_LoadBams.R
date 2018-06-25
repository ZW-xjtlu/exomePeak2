library(testthat)

test_that( "Check Bams", {
  MeRIP_Seq_Alignment <- scan_merip_bams(bam_ip = c("./bam/SRR1182619.bam",
                                             "./bam/SRR1182621.bam",
                                             "./bam/SRR1182623.bam"),
                              bam_input = c("./bam/SRR1182620.bam",
                                            "./bam/SRR1182622.bam",
                                            "./bam/SRR1182624.bam"),
                              bam_treated_ip = c("./bam/SRR1182603.bam",
                                                 "./bam/SRR1182605.bam"),
                              bam_treated_input = c("./bam/SRR1182604.bam",
                                                    "./bam/SRR1182606.bam"),
                              paired_end = T,
                              bam_files = NULL,
                              design_ip = NULL,
                              design_treatment = NULL)

  MeRIP_Seq_Alignment2 <- scan_merip_bams(bam_files =
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
                                  design_ip = metadata(MeRIP_Seq_Alignment[[1]])$design_IP,
                                  design_treatment = metadata(MeRIP_Seq_Alignment[[1]])$design_Treatment,
                                  paired_end = T)

  expect_that(MeRIP_Seq_Alignment, is_a("list") ) #Later when we use OOP, we should change the list into our own S4 object.
  expect_that(MeRIP_Seq_Alignment, is_identical_to(MeRIP_Seq_Alignment) )
} )


