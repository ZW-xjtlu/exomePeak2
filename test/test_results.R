library(testthat)

test_that( "test results export of exomePeak2", {

SEP_meth <- readRDS("SEP_meth.rds")
SEP_meth <- readRDS("./Developing/SEP_meth.rds")
results_df <- Results(SEP_meth)
exportResults(SEP_meth)

})
