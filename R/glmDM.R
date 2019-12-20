#' @title Estimation and Inference on RNA Differential Modification Log2 Fold Changes with Generalized Linear Model.
#'
#' @param sep a \code{\link{SummarizedExomePeak}} object.
#'
#' @param glm_type a \code{character} speciefies the type of Generalized Linear Model (GLM) fitted for the purpose of statistical inference during peak calling, which can be one of the \code{c("DESeq2", "NB", "Poisson")}.
#'
#' \describe{
#' \item{\strong{\code{DESeq2}}}{Fit the GLM defined in the function \code{\link{DESeq}}, which is the NB GLM with regulated estimation of the overdispersion parameters.}
#'
#' \item{\strong{\code{NB}}}{Fit the Negative Binomial (NB) GLM.}
#'
#' \item{\strong{\code{Poisson}}}{Fit the Poisson GLM.}
#' }
#'
#' By default, the DESeq2 GLMs are fitted on the data set with > 1 biological replicates for both the IP and input samples, the Poisson GLM will be fitted otherwise.
#'
#' @param LFC_shrinkage a \code{character} for the method of emperical bayes shrinkage on log2FC, could be one of \code{c("apeglm", "ashr", "none")}; Default \code{= "apeglm"}.
#'
#' see \code{\link{lfcShrink}} for details; if "none" is selected, only the MLE will be returned.
#'
#' @param ... Optional arguments passed to \code{\link{DESeq}}
#'
#' @description \code{glmDM} perform inference and estimation on RNA differential modification log2FC.
#'
#' GLMs with interactive design between dummy variables of IP/input and Treatment/control are fitted for each peaks/sites:
#'
#' \deqn{log2(Q) = intercept + I(Treatment) + I(IP) + I(IP)*I(Treatment)}
#'
#' The log2FC and the associated statistics are based on the coefficient estimate of the interactive term: \eqn{I(IP)*I(Treated)}.
#'
#' Under default setting, the returned log2FC are the RR estimates with Couchey priors defined in \code{\link{apeglm}}.
#'
#' @import SummarizedExperiment
#' @import apeglm
#' @import DESeq2
#'
#' @docType methods
#'
#' @return a \code{SummarizedExomPeak} object.
#'
#' @examples
#'
#' ### Load the example SummarizedExomPeak object
#' f1 = system.file("extdata", "sep_ex_dm.rds", package="exomePeak2")
#'
#' sep <- readRDS(f1)
#'
#' ### Normalize the GC contents biases
#' sep <- normalizeGC(sep)
#'
#' ### Calculate GLM Statistics on the Modification Peaks
#' sep <- glmDM(sep)
#'
#' @aliases glmDM
#'
#' @rdname glmDM-methods
#'
#' @seealso \code{\link{glmM}}
#'
#' @export

setMethod("glmDM",
          "SummarizedExomePeak",
           function(sep,
                    glm_type = c("auto","Poisson", "NB", "DESeq2"),
                    LFC_shrinkage = c("apeglm","ashr","none"),
                    ...) {

  LFC_shrinkage = match.arg(LFC_shrinkage)

  glm_type = match.arg(glm_type)

   if( !(any(sep$design_Treatment) & any(!sep$design_Treatment) ) ){
     stop("Your data has no interactive design / treatment groups, please use glm_M().\n")
   }

  if(glm_type == "auto") {
    if( all( table(colData(sep)$design_IP,colData(sep)$design_Treatment) > 1 ) ) {
      glm_type <- "DESeq2"
    } else {
      glm_type <- "Poisson"
    }
  }

  if(glm_type == "Poisson") {
    message("Differential modification analysis with Poisson GLM ... ", appendLF = FALSE)
  }

  if(glm_type == "NB") {
    message("Differential modification analysis with NB GLM ... ", appendLF = FALSE)
  }

  if(glm_type == "DESeq2") {
    message("Differential modification analysis with regularized NB GLM ... ", appendLF = FALSE)
  }

  if(is.null(colData( sep )$sizeFactor)) {

    sep <- estimateSeqDepth(sep)

  }

  indx_mod <- grepl("peak_", rownames( sep ) )

  SE_M <- sep

  SE_M$IPinput = "input"

  SE_M$IPinput[SE_M$design_IP] = "IP"

  SE_M$IPinput = factor(SE_M$IPinput)

  SE_M$Perturbation = "Control"

  SE_M$Perturbation[SE_M$design_Treatment] = "Treatment"

  SE_M$Perturbation = factor( SE_M$Perturbation )

  if(!is.null(GCsizeFactors( sep ))) {

    gc_na_indx <- rowSums( is.na(GCsizeFactors(sep)) ) > 0

    Cov = ~ Perturbation + IPinput + Perturbation:IPinput

    dds = suppressMessages( DESeqDataSet(se = SE_M[(!gc_na_indx) & indx_mod,], design = Cov) )

    glm_off_sets <- GCsizeFactors(sep)[(!gc_na_indx) & indx_mod,]

    #Normalization to make the row geometric means = 0 (since DESeq2 only cares about the difference)
    #and this norm factor is still under the original scale (not log scale glm off set).

    centered_off_sets <- exp(glm_off_sets) / exp(rowMeans(glm_off_sets))

    normalizationFactors(dds) <- centered_off_sets

    rm(glm_off_sets,centered_off_sets)

  } else {

    Cov = ~ Perturbation + IPinput + Perturbation:IPinput

    dds = suppressMessages( DESeqDataSet(se = SE_M[indx_mod,], design = Cov) )

  }

  dds$IPinput <- relevel(dds$IPinput, "input")

  dds$Perturbation <- relevel(dds$Perturbation, "Control")

  if(glm_type == "Poisson"){
    dispersions(dds) = 0
  }

  if(glm_type == "NB"){
    dds =  suppressMessages( estimateDispersions( dds, fitType = "mean" ) )
  }

  if(glm_type == "DESeq2"){
    dds = suppressMessages( estimateDispersions( dds ) )
  }


  dds = suppressMessages( nbinomWaldTest( dds ) )

  #Generation of the DESeq2 report.
  DS_result <- as.data.frame( suppressMessages( results( dds ) ) )

  message("OK")

  DS_result <- DS_result[,c("log2FoldChange","lfcSE","pvalue","padj")]
  colnames(DS_result) <- c("log2fcDiffMod.MLE","log2fcDiffMod.MLE.SE","log2fcDiffMod.pvalue","log2fcDiffMod.padj")

  #Include reads count
  DS_result$ReadsCount.IP.Treated <- rowSums( cbind(assay(dds)[,colData(dds)$design_IP & colData(dds)$design_Treatment]) )
  DS_result$ReadsCount.input.Treated <- rowSums( cbind(assay(dds)[,!colData(dds)$design_IP & colData(dds)$design_Treatment] ))
  DS_result$ReadsCount.IP.Control <- rowSums( cbind(assay(dds)[,colData(dds)$design_IP & !colData(dds)$design_Treatment] ))
  DS_result$ReadsCount.input.Control <- rowSums( cbind(assay(dds)[,!colData(dds)$design_IP & !colData(dds)$design_Treatment] ))

  #Calculate estimates of other contrasts
  Expr_Control_design_MLE <- as.data.frame( suppressWarnings( suppressMessages( results( dds, contrast = c(1,0,0,0)) ) ))
  DS_result$log2Expr.Control.MLE <- Expr_Control_design_MLE[,"log2FoldChange"]
  DS_result$log2Expr.Control.MLE.SE <- Expr_Control_design_MLE[,"lfcSE"]
  rm(Expr_Control_design_MLE)

  Expr_Treated_design_MLE <- as.data.frame( suppressWarnings( suppressMessages( results( dds, contrast = c(1,1,0,0)) ) ))
  DS_result$log2Expr.Treated.MLE <- Expr_Treated_design_MLE[,"log2FoldChange"]
  DS_result$log2Expr.Treated.MLE.SE <- Expr_Treated_design_MLE[,"lfcSE"]
  rm(Expr_Treated_design_MLE)

  DiffExpr_design_MLE <- as.data.frame( suppressWarnings( suppressMessages( results( dds, contrast = c(0,1,0,0)) ) ))
  DS_result$log2fcDiffExpr.MLE <- DiffExpr_design_MLE[,"log2FoldChange"]
  DS_result$log2fcDiffExpr.MLE.SE <- DiffExpr_design_MLE[,"lfcSE"]
  DS_result$log2fcDiffExpr.pvalue <- DiffExpr_design_MLE[,"pvalue"]
  DS_result$log2fcDiffExpr.padj <- DiffExpr_design_MLE[,"padj"]
  rm(DiffExpr_design_MLE)

  Mod_Control_design_MLE <- as.data.frame( suppressWarnings( suppressMessages( results( dds, contrast = c(0,0,1,0)) ) ))
  DS_result$log2fcMod.Control.MLE <- Mod_Control_design_MLE[,"log2FoldChange"]
  DS_result$log2fcMod.Control.MLE.SE  <- Mod_Control_design_MLE[,"lfcSE"]
  DS_result$log2fcMod.Control.pvalue <- Mod_Control_design_MLE[,"pvalue"]
  DS_result$log2fcMod.Control.padj <- Mod_Control_design_MLE[,"padj"]
  rm(Mod_Control_design_MLE)

  Mod_Treated_design_MLE <- as.data.frame( suppressWarnings( suppressMessages( results( dds, contrast = c(0,0,1,1)) ) ))
  DS_result$log2fcMod.Treated.MLE <- Mod_Treated_design_MLE[,"log2FoldChange"]
  DS_result$log2fcMod.Treated.MLE.SE  <- Mod_Treated_design_MLE[,"lfcSE"]
  DS_result$log2fcMod.Treated.pvalue <- Mod_Treated_design_MLE[,"pvalue"]
  DS_result$log2fcMod.Treated.padj <- Mod_Treated_design_MLE[,"padj"]
  rm(Mod_Treated_design_MLE)

  AbsDiffMod_design_MLE <- as.data.frame( suppressWarnings( suppressMessages( results( dds, contrast = c(0,1,0,1)) ) ))
  DS_result$log2fcAbsDiffMod.MLE <- AbsDiffMod_design_MLE[,"log2FoldChange"]
  DS_result$log2fcAbsDiffMod.MLE.SE  <- AbsDiffMod_design_MLE[,"lfcSE"]
  DS_result$log2fcAbsDiffMod.pvalue <- AbsDiffMod_design_MLE[,"pvalue"]
  DS_result$log2fcAbsDiffMod.padj <- AbsDiffMod_design_MLE[,"padj"]
  rm(AbsDiffMod_design_MLE)

  #Calculate additional MAP estimates if LFCs are set != "none"
  if(LFC_shrinkage != "none" & glm_type != "Poisson") {

    DS_result$log2fcDiffExpr.MAP <- as.data.frame( suppressMessages( lfcShrink( dds=dds, coef = 2, type = LFC_shrinkage  ) ) )$log2FoldChange
    DS_result$log2fcDiffMod.MAP <- as.data.frame( suppressMessages( lfcShrink( dds=dds, coef = 4, type = LFC_shrinkage  ) ) )$log2FoldChange
    DS_result$log2FoldChange <- DS_result$log2fcDiffMod.MAP

    DS_result <- DS_result[,c("ReadsCount.input.Control","ReadsCount.IP.Control",
                              "ReadsCount.input.Treated","ReadsCount.IP.Treated",
                              "log2Expr.Control.MLE","log2Expr.Control.MLE.SE",
                              "log2Expr.Treated.MLE","log2Expr.Treated.MLE.SE",
                              "log2fcMod.Control.MLE","log2fcMod.Control.MLE.SE",
                              "log2fcMod.Control.pvalue","log2fcMod.Control.padj",
                              "log2fcMod.Treated.MLE","log2fcMod.Treated.MLE.SE",
                              "log2fcMod.Treated.pvalue","log2fcMod.Treated.padj",
                              "log2fcDiffExpr.MLE","log2fcDiffExpr.MLE.SE","log2fcDiffExpr.MAP",
                              "log2fcDiffExpr.pvalue","log2fcDiffExpr.padj",
                              "log2fcAbsDiffMod.MLE","log2fcAbsDiffMod.MLE.SE",
                              "log2fcAbsDiffMod.pvalue","log2fcAbsDiffMod.padj",
                              "log2fcDiffMod.MLE","log2fcDiffMod.MLE.SE","log2fcDiffMod.MAP",
                              "log2fcDiffMod.pvalue","log2fcDiffMod.padj")]

    major_statistics <- DS_result[,c("log2fcDiffMod.MAP","log2fcDiffMod.pvalue","log2fcDiffMod.padj")]

    colnames(major_statistics) = c("log2FoldChange","pvalue","padj")

    DS_result <- cbind(DS_result, major_statistics)

    rm(major_statistics)

  } else {
    DS_result <- DS_result[,c("ReadsCount.input.Control","ReadsCount.IP.Control",
                              "ReadsCount.input.Treated","ReadsCount.IP.Treated",
                              "log2Expr.Control.MLE","log2Expr.Control.MLE.SE",
                              "log2Expr.Treated.MLE","log2Expr.Treated.MLE.SE",
                              "log2fcMod.Control.MLE","log2fcMod.Control.MLE.SE",
                              "log2fcMod.Control.pvalue","log2fcMod.Control.padj",
                              "log2fcMod.Treated.MLE","log2fcMod.Treated.MLE.SE",
                              "log2fcMod.Treated.pvalue","log2fcMod.Treated.padj",
                              "log2fcDiffExpr.MLE","log2fcDiffExpr.MLE.SE",
                              "log2fcDiffExpr.pvalue","log2fcDiffExpr.padj",
                              "log2fcAbsDiffMod.MLE","log2fcAbsDiffMod.MLE.SE",
                              "log2fcAbsDiffMod.pvalue","log2fcAbsDiffMod.padj",
                              "log2fcDiffMod.MLE","log2fcDiffMod.MLE.SE",
                              "log2fcDiffMod.pvalue","log2fcDiffMod.padj")]

    major_statistics <- DS_result[,c("log2fcDiffMod.MLE","log2fcDiffMod.pvalue","log2fcDiffMod.padj")]

    colnames(major_statistics) = c("log2FoldChange","pvalue","padj")

    DS_result <- cbind(DS_result, major_statistics)

    rm(major_statistics)

  }

  if (!is.null(GCsizeFactors( sep ))) {

     DS_final_rst <- matrix( NA, nrow = nrow(SE_M[indx_mod,]), ncol = ncol(DS_result) )

    colnames( DS_final_rst) <- colnames(DS_result)

     DS_final_rst <- as.data.frame( DS_final_rst)

     DS_final_rst[(!gc_na_indx)[indx_mod],] <- as.data.frame( DS_result )

     rm(DS_result)

  } else {

    DS_final_rst <- DS_result

    rm(DS_result)

  }

  rownames(  DS_final_rst ) = rownames( SE_M )[indx_mod]

  exomePeak2Results( sep ) = as.data.frame(  DS_final_rst )

  return( sep )

})
