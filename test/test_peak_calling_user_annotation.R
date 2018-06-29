library(testthat)

test_that( "Peak Calling with user provided annotation", {
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
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

  hg19_miCLIP <- readRDS("hg19_miCLIP.rds")

  SummarizedExomePeaks <- merip_peak_calling( MeRIP_Seq_Alignment,
                                              txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                              provided_annotation = hg19_miCLIP )

  SummarizedExomePeaks <- estimate_size_factors(SummarizedExomePeaks)

  SummarizedExomePeaks <- GC_normalization(SummarizedExomePeaks, bsgenome = Hsapiens)

  SummarizedExomePeaks <- glm_dm(SummarizedExomePeaks)

  plot_size_factors(SummarizedExomePeaks)

  guitar_plot(SummarizedExomePeaks,
              save_pdf_prefix = "DM_SB",
              guitar_coordinate = readRDS("/Users/zhenwei/Datasets/Gtcoords/Gtcoord_hg19.rds"))

  exon_length_plot(SummarizedExomePeaks,
                   txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                   save_pdf_prefix = "DM_SB")

  aabundance_GC_plot(SummarizedExomePeaks,
                    bsgenome = Hsapiens,
                    save_pdf_prefix = "DM_SB")

  beta_GC_plot(
    SummarizedExomePeaks,
    bsgenome = Hsapiens,
    save_pdf_prefix = "Meth_SB"
  )

  abundance_GC_plot(SummarizedExomePeaks,
                    bsgenome = Hsapiens,
                    save_pdf_prefix = "Meth_SB")

  saveRDS(SummarizedExomePeaks, "SEP_sb.rds")


} )
