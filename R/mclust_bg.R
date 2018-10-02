#'@title Find background of merip-seq signal with model based clustering.
#'@param se_peak_counts A SummarizedExperiment object.
#'@param alpha The small offset added when calculating the M value.
#'@description The function generate background index on the rowData of the summarizedExperiment.
#'@return A dummy variable of the background index.
#'@details This function will do the following jobs:
#'
#'1. Filter the rows (modification sites) by average count.
#'
#'2. Fit multivariate gaussian mixture model with 2 mixing component to differentiate background and biological signal.
#'
#'- depend on whether the bsgenome is provided, the values GC content will be used as one of the dimention of the clustering.
#'
#'3. Classify the background and signal using bayes classifier by the learned model.
#'
#'4. Return the index for the bins that are classified into background.
#'
#'@import BSgenome
#'@import SummarizedExperiment
#'@import mclust
#'
mclust_bg <- function(se_peak_counts,
                      alpha = 1){


#1. Calculate the methylation levels as M values.
mean_input_control <- rowMeans(assay(se_peak_counts)[rowData(se_peak_counts)$indx_gc_est,
                                                     !(se_peak_counts$design_IP)&!(se_peak_counts$design_Treatment)]) + alpha

mean_input_treated <- rowMeans(assay(se_peak_counts)[rowData(se_peak_counts)$indx_gc_est,
                                                     !(se_peak_counts$design_IP)&(se_peak_counts$design_Treatment)]) + alpha

IP_control <- assay(se_peak_counts)[rowData(se_peak_counts)$indx_gc_est,
                                    (se_peak_counts$design_IP)&!(se_peak_counts$design_Treatment)] + alpha

IP_treated <- assay(se_peak_counts)[rowData(se_peak_counts)$indx_gc_est,
                                    (se_peak_counts$design_IP)&(se_peak_counts$design_Treatment)] + alpha

M_value_control <- IP_control/mean_input_control
M_value_treated <- IP_treated/mean_input_treated

rm(mean_input_control, mean_input_treated, IP_control, IP_treated)

#2. Apply multivariate gaussian mixture model on all M levels and the GC contents.
model_matrix <- log(cbind(M_value_control,M_value_treated))
rm(M_value_control,M_value_treated)

colnames(model_matrix) <- paste0("sample_",1:ncol(model_matrix))
model_matrix <- as.data.frame(model_matrix)
if(!is.null(rowData(se_peak_counts)$gc_contents)) model_matrix$GC <- rowData(se_peak_counts)$gc_contents[rowData(se_peak_counts)$indx_gc_est]
mod_mix <- densityMclust(model_matrix, G = 2, modelNames = "VVV")

#3. Classify the bins using bayesian classifier
bg_class <- which.min( colSums( mod_mix$parameters$mean[seq_len(sum(grepl("sample",colnames(model_matrix) ))),] ) )

indx_bg <- vector(length = nrow(se_peak_counts), mode = "logical")

indx_bg[rowData(se_peak_counts)$indx_gc_est] <- mod_mix$classification == bg_class

rm(bg_class)

return(indx_bg)
}


