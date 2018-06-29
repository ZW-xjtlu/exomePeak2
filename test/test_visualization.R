library(testthat)

test_that( "Visualization", {
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

  SEP_dm <- readRDS("SEP_dm2.rds")
  SEP_sb <- readRDS("SEP_sb.rds")

  plot_size_factors(SEP_dm)
  plot_size_factors(SEP_sb)

  SEP_dm <- estimate_size_factors(SEP_dm)
  SEP_meth <- estimate_size_factors(SEP_meth)
  SEP_dm <- GC_normalization(SEP_dm, bsgenome = Hsapiens)
  SEP_meth <- GC_normalization(SEP_meth,bsgenome = Hsapiens)
  SEP_meth <- glm_meth(SEP_meth)
  SEP_dm <- glm_dm(SEP_dm)

  plot_size_factors(SEP_meth)
  plot_size_factors(SEP_dm)

  guitar_plot(SEP_meth,
              save_pdf_prefix = "Meth_mcqn",
              guitar_coordinate = readRDS("/Users/zhenwei/Datasets/Gtcoords/Gtcoord_hg19.rds"))


  guitar_plot(SEP_dm,
              save_pdf_prefix = "DM_mcqn",
              guitar_coordinate = readRDS("/Users/zhenwei/Datasets/Gtcoords/Gtcoord_hg19.rds"))


  exon_length_plot(SEP_meth,
                   txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                   save_pdf_prefix = "Meth_mcqn")

  exon_length_plot(SEP_dm,
                   txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                   save_pdf_prefix = "DM_mcqn")

  abundance_GC_plot(SEP_meth,
                bsgenome = Hsapiens,
                save_pdf_prefix = "Meth_mcqn")

  abundance_GC_plot(SEP_dm,
                bsgenome = Hsapiens,
                save_pdf_prefix = "DM_mcqn")

  SEP_meth <- GC_normalization(SEP_meth,
                               bsgenome = Hsapiens,
                               scope = "methylation")

  SEP_dm <- GC_normalization(SEP_dm,
                                bsgenome = Hsapiens,
                               scope = "methylation")

  abundance_GC_plot(SEP_meth,
                    bsgenome = Hsapiens,
                    save_pdf_prefix = "Meth_mcqn")

  abundance_GC_plot(SEP_dm,
                    bsgenome = Hsapiens,
                    save_pdf_prefix = "DM_mcqn")

  beta_GC_plot(
    SEP_meth,
    bsgenome = Hsapiens,
    save_pdf_prefix = "Meth_mcqn"
  )

  beta_GC_plot(
    SEP_dm,
    bsgenome = Hsapiens,
    save_pdf_prefix = "DM_mcqn"
  )

  SEP_dm2 <- readRDS( "SEP_dm2.rds" )

  SEP_dm2 <- estimate_size_factors(SEP_dm2)

  plot_size_factors(SEP_dm2)

  SEP_dm2 <- GC_normalization(SEP_dm2, bsgenome = Hsapiens)

  SEP_dm2 <- glm_dm(SEP_dm2)

  beta_GC_plot(
  SEP_dm2,
  bsgenome = Hsapiens,
  save_pdf_prefix = "DM_2"
  )

  exon_length_plot(
    SEP_dm2,
    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
    save_pdf_prefix = "DM_2"
  )

  guitar_plot(SEP_dm2,
              save_pdf_prefix = "DM2",
              guitar_coordinate = readRDS("/Users/zhenwei/Datasets/Gtcoords/Gtcoord_hg19.rds"))

  abundance_GC_plot(SEP_dm2,
                    bsgenome = Hsapiens,
                    save_pdf_prefix = "DM2")
} )

