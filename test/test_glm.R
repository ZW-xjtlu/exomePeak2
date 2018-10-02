library(testthat)

test_that( "Quantification", {
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(exomePeak2)

  SEP_meth <- readRDS("SEP_meth.rds")
  SEP_meth <- estimateSeqDepth(SEP_meth)
  plotSizeFactors(SEP_meth)
  SEP_meth <- normalizGC(SEP_meth,bsgenome = Hsapiens)
  SEP_meth <- glmMeth(SEP_meth)
  SEP_meth <- glmDM(SEP_meth)

  saveRDS(SEP_meth,"SEP_meth.rds")

} )
