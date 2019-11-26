#'@title Find background of merip-seq signal with model based clustering.
#'@param se_peak_counts A SummarizedExperiment object.
#'@param alpha The small offset added when calculating the M value.
#'@description The function generate background index on the rowData of the summarizedExperiment.
#'@return A dummy variable of the background index.
#'@details This function will do the following jobs:
#'
#'1. Filter the rows (modification sites) by average count.
#'
#'2. Fit multivariate Gaussian mixture model with 2 mixing component to differentiate background and biological signal.
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
#'
mclust_bg <- function(se_peak_counts,
                      alpha = 1){

#1. Calculate the methylation levels as M values.

indx_input_control <- !(se_peak_counts$design_IP) & !(se_peak_counts$design_Treatment)

if(sum(indx_input_control) == 1){
mean_input_control <- assay(se_peak_counts)[rowData(se_peak_counts)$indx_gc_est,indx_input_control] + alpha
}else{
mean_input_control <- rowMeans(assay(se_peak_counts)[rowData(se_peak_counts)$indx_gc_est,indx_input_control]) + alpha
}

rm(indx_input_control)

indx_input_treated <- !(se_peak_counts$design_IP) & (se_peak_counts$design_Treatment)

if(!any(indx_input_treated)){
}else{
  if(sum(indx_input_treated) == 1){
    mean_input_treated <- assay(se_peak_counts)[rowData(se_peak_counts)$indx_gc_est,indx_input_treated] + alpha
  } else {
    mean_input_treated <- rowMeans(assay(se_peak_counts)[rowData(se_peak_counts)$indx_gc_est,indx_input_treated]) + alpha
  }
}

IP_control <- assay(se_peak_counts)[rowData(se_peak_counts)$indx_gc_est,
                                    (se_peak_counts$design_IP)&!(se_peak_counts$design_Treatment)] + alpha

indx_IP_treated <- (se_peak_counts$design_IP)&(se_peak_counts$design_Treatment)

if(any(indx_IP_treated)){
 IP_treated <- assay(se_peak_counts)[rowData(se_peak_counts)$indx_gc_est,indx_IP_treated] + alpha
}

M_value_control <- IP_control/mean_input_control

if(any(c(indx_input_treated, indx_IP_treated))) {
M_value_treated <- IP_treated/mean_input_treated
model_matrix <- log(cbind(M_value_control,M_value_treated))
rm(M_value_control,M_value_treated, mean_input_treated, IP_treated)
} else{
model_matrix <- log(M_value_control)
rm(M_value_control)
}

rm(mean_input_control, IP_control, indx_IP_treated, indx_input_treated)

#2. Apply multivariate Gaussian mixture model on all M levels and the GC contents.

model_matrix <- as.data.frame(model_matrix)
colnames(model_matrix) <- paste0("sample_",seq_len(ncol(model_matrix)))
if(!is.null(rowData(se_peak_counts)$gc_contents)) model_matrix$GC <- rowData(se_peak_counts)$gc_contents[rowData(se_peak_counts)$indx_gc_est]

#Apply univariate Gaussian mixture if there is only one column.
#Classify the bins using bayesian classifier
if(ncol(model_matrix) == 1){

mod_mix <- densityMclust(as.numeric( model_matrix[,1]), G = 2, modelNames = "V")
bg_class <- which.min( mod_mix$parameters$mean )

} else {

mod_mix <- densityMclust(model_matrix, G = 2, modelNames = "VVV")
par_samples <- mod_mix$parameters$mean[seq_len(sum(grepl("sample",colnames(model_matrix) ))),]

if(ncol(model_matrix) == 2){
  bg_class <- which.min( par_samples )
}else{
  bg_class <- which.min( colSums( par_samples  ) )
}

rm(par_samples)

}

indx_bg <- vector(length = nrow(se_peak_counts), mode = "logical")

indx_bg[rowData(se_peak_counts)$indx_gc_est] <- mod_mix$classification == bg_class

rm(bg_class)

return(indx_bg)
}


