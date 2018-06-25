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

  SEP_dm <- readRDS("SEP_dm.rds")

  SEP_meth <- readRDS("SEP_meth.rds")

  SEP_dm <- estimate_size_factors(SEP_dm)
  SEP_meth <- estimate_size_factors(SEP_meth)

  plot_size_factors(SEP_meth)
  plot_size_factors(SEP_dm)

  SEP_dm <- GC_normalization(SEP_dm, bsgenome = Hsapiens)
  SEP_meth <- GC_normalization(SEP_meth,bsgenome = Hsapiens)

  SEP_meth_CQN <- glm_meth(SEP_meth)

  SEP_meth <- readRDS("SEP_meth.rds")

  SEP_meth_noCQN <- glm_meth(SEP_meth)

  sum(SEP_meth_CQN$DESeq2Result$padj < .05)
  sum(SEP_meth_noCQN$DESeq2Result$padj < .05)

  dir_gtcoord = "/Users/zhenwei/Datasets/Gtcoords/Gtcoord_hg19.rds"
  require(Guitar)
  GuitarPlot(list(CQN =  rowRanges(SEP_meth_CQN$SE)[grepl("meth",rownames(SEP_meth_CQN$SE))][SEP_meth_CQN$DESeq2Result$padj < .05],
                  NO_CQN = rowRanges(SEP_meth_noCQN$SE)[grepl("meth",rownames(SEP_meth_CQN$SE))][SEP_meth_noCQN$DESeq2Result$padj < .05]),
             readRDS(dir_gtcoord))


  SEP_dm_CQN <- glm_dm(SEP_dm)

  length(which( SEP_dm_CQN$DESeq2Result$padj < 0.05) )

  length(which( SEP_dm_CQN$DESeq2Result$pvalue < 0.05 &
                  SEP_dm_CQN$DESeq2Result$log2FoldChange < 0) )


  SEP_dm <- readRDS("SEP_dm.rds")
  SEP_dm_noCQN <- glm_dm(SEP_dm)

  length(which( SEP_dm_noCQN$DESeq2Result$padj < 0.05) )

  length(which( SEP_dm_noCQN$DESeq2Result$pvalue < 0.05 &
                  SEP_dm_noCQN$DESeq2Result$log2FoldChange < 0) )

  #One thing worth to know that, using which shrinkage on estimate does not change p values (wald test),
  #which still use the not shrinked estimate (not ridge regression p value, but regular wald test).


  GuitarPlot(list(CQN =  rowRanges(SEP_dm_CQN$SE)[grepl("meth",rownames(SEP_dm_CQN$SE))][which( SEP_dm_CQN$DESeq2Result$pvalue < .05 &
                                                                                                  SEP_dm_CQN$DESeq2Result$log2FoldChange < 0)],
                  NO_CQN = rowRanges(SEP_dm_noCQN$SE)[grepl("meth",rownames(SEP_dm_CQN$SE))][which( SEP_dm_noCQN$DESeq2Result$pvalue < .05 &
                                                                                                      SEP_dm_noCQN$DESeq2Result$log2FoldChange < 0)]),
             readRDS(dir_gtcoord))


  p_threshold_CQN <- meripQC::Decision_infresult(SEP_dm_CQN$DESeq2Result,
                              Min_mod = 3000,
                              P_cut = 0.01,
                              DM = TRUE,
                              Exp_dir = "hypo",
                              HDER = "with CQN")$Cut_Val_expected

  p_threshold_noCQN <- meripQC::Decision_infresult(SEP_dm_noCQN$DESeq2Result,
                                                 Min_mod = 3000,
                                                 P_cut = 0.01,
                                                 DM = TRUE,
                                                 Exp_dir = "hypo",
                                                 HDER = "with CQN")$Cut_Val_expected

  GuitarPlot(list(CQN =  rowRanges(SEP_dm_CQN$SE)[grepl("meth",rownames(SEP_dm_CQN$SE))][which( SEP_dm_CQN$DESeq2Result$pvalue < p_threshold_CQN &
                                                                                                  SEP_dm_CQN$DESeq2Result$log2FoldChange < 0)],
                  NO_CQN = rowRanges(SEP_dm_noCQN$SE)[grepl("meth",rownames(SEP_dm_CQN$SE))][which( SEP_dm_noCQN$DESeq2Result$pvalue < p_threshold_noCQN &
                                                                                                      SEP_dm_noCQN$DESeq2Result$log2FoldChange < 0)]),
             readRDS(dir_gtcoord))


  expect_that( SummarizedExomePeaks, is_a("list") ) #Later when we use OOP, we should change the list into our own S4 object.
} )
