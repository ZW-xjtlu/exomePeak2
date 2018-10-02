library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(exomePeak2)
library(BSgenome.Hsapiens.UCSC.hg19)

SEP_dm2 <- readRDS("SEP_dm2.rds")
SEP_sb <- readRDS("SEP_sb.rds")
SEP_dm2$SE$sizeFactor <- SEP_sb$SE$sizeFactor #Use size factor estimated from single based data.

SEP_dm2 <- GC_normalization(SEP_dm2)

SEP_dm2 <- glm_dm(SEP_dm2)

plot_size_factors(SEP_dm2)

guitar_plot(SEP_dm2,
            save_pdf_prefix = "DM3",
            guitar_coordinate = readRDS("/Users/zhenwei/Datasets/Gtcoords/Gtcoord_hg19.rds"))

exon_length_plot(SEP_dm2,
                 txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                 save_pdf_prefix = "DM3")

aabundance_GC_plot(SEP_dm2,
                   bsgenome = Hsapiens,
                   save_pdf_prefix = "DM3")

beta_GC_plot(
  SEP_dm2,
  bsgenome = Hsapiens,
  save_pdf_prefix = "DM3"
)

abundance_GC_plot(SEP_dm2,
                  bsgenome = Hsapiens,
                  save_pdf_prefix = "DM3")

#Conclusion: the size factor selection method is really not optimal.
#The peak refinement is good, but the peak calling is shit,

SEP_dm2$SE$sizeFactor <- size_factors_cdf

size_factors_cdf

size_factor_M <- estimate_size_factors( SEP_dm2 ,"Methylation")$SE$sizeFactor

ratios <- size_factor_M/size_factors_cdf

ratios[c(1,2,3)] #IP enrichment

ratios[c(7,8)] #Treatment enrichment



SEP_dm2 <- GC_normalization(SEP_dm2)

SEP_dm2 <- glm_dm(SEP_dm2)

plot_size_factors(SEP_dm2)

guitar_plot(SEP_dm2,
            save_pdf_prefix = "DM4",
            guitar_coordinate = readRDS("/Users/zhenwei/Datasets/Gtcoords/Gtcoord_hg19.rds"))

exon_length_plot(SEP_dm2,
                 txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                 save_pdf_prefix = "DM4")

aabundance_GC_plot(SEP_dm2,
                   bsgenome = Hsapiens,
                   save_pdf_prefix = "DM4")

beta_GC_plot(
  SEP_dm2,
  bsgenome = Hsapiens,
  save_pdf_prefix = "DM4"
)

abundance_GC_plot(SEP_dm2,
                  bsgenome = Hsapiens,
                  save_pdf_prefix = "DM4")

#Conclusion: this makes no difference from our background method....


