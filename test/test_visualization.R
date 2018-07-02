library(testthat)

test_that( "Visualization", {
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(exomePeak2)

  SEP_meth <- readRDS("SEP_meth.rds")

  plotSizeFactors(SEP_meth)

  plotGuitar(SEP_meth,
              save_pdf_prefix = "Meth_mcqn",
              guitar_coordinate = readRDS("/Users/zhenwei/Datasets/Gtcoords/Gtcoord_hg19.rds"))

  plotExonLength(SEP_meth,
                   txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                   save_pdf_prefix = "Meth_mcqn")


  plotReadsGC(
    SEP_meth,
    bsgenome = Hsapiens,
    save_pdf_prefix = "Meth_mcqn"
  )

  plotBetaGC(
    SEP_meth,
    bsgenome = Hsapiens,
    save_pdf_prefix = "Meth_mcqn"
  )

} )

