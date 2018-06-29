library(testthat)

test_that( "test results export of exomePeak2", {

SEP_sb <- readRDS("SEP_sb.rds")
SEP_meth <- readRDS("./Developing/SEP_meth.rds")

save_results(
              SEP_sb,
              format = "txt",
              cut_off_padj = 0.05,
              cut_off_log2FC = 1,
              min_num_of_positive = 2000,
              inhibit_filter = TRUE
)

save_results(
             SEP_sb,
             format = "RDS",
             cut_off_padj = 0.01,
             cut_off_log2FC = 1,
             min_num_of_positive = 2000,
             expected_direction = "hypo"
)

SEP_meth <- glm_meth(SEP_meth)

save_results(
             SEP_meth,
             format = "txt",
             cut_off_padj = 0.05,
             cut_off_log2FC = 1,
             min_num_of_positive = 2000,
             inhibit_filter = TRUE
)

save_results(
             SEP_meth,
             format = "RDS",
             cut_off_padj = 0.01,
             cut_off_log2FC = 1,
             min_num_of_positive = 2000,
             expected_direction = "hypo"
)

save_results(
             SEP_meth,
             format = "BED",
             cut_off_padj = 0.01,
             cut_off_log2FC = 1,
             min_num_of_positive = 2000,
             expected_direction = "hypo"
)


save_results(
             SEP_meth,
             format = "txt",
             cut_off_padj = 0.001,
             cut_off_log2FC = 1,
             min_num_of_positive = 2000,
             expected_direction = "hypo"
)


})
