#' @title Estimation and Inference on RNA Differential Modification Log2 Fold Changes with Generalized Linear Model.
#'
#' @param sep a \code{\link{summarizedExomePeak}} object.
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
#' see \code{\link{lfcShrink}} for more details.
#'
#' @param ... Optional arguments passed to \code{\link{DESeq}}
#'
#' @description \code{glmDM} perform inference and estimation on RNA differential modification log2FC.
#'
#' GLMs with interactive design between dummy variables of IP/input and Treatment/control are fitted for each peaks/sites:
#'
#' \deqn{log2(Q) = intercept + I(Treatment) + I(IP) + I(IP)*I(Treatment)}
#'
#' The log2FC and the associated statistics are based on the coefficient estimate of the interactive term： \eqn{I(IP)*I(Treated)}.
#'
#' Under default setting, the returned log2FC are the RR estimates with Couchey priors defined in \code{\link{apeglm}}.
#'
#' @import SummarizedExperiment
#'
#' @import DESeq2
#'
#' @docType methods
#'
#' @examples
#'
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#'
#' aln <- scanMeripBAM(
#' bam_ip = c("IP_rep1.bam",
#'            "IP_rep2.bam",
#'            "IP_rep3.bam"),
#' bam_input = c("input_rep1.bam",
#'               "input_rep2.bam",
#'               "input_rep3.bam"),
#' bam_treated_ip = c("IP_treated_rep1.bam",
#'                    "IP_treated_rep2.bam"),
#' bam_treated_input = c("input_treated_rep1.bam",
#'                       "input_treated_rep2.bam"),
#' paired_end = TRUE
#' )
#'
#' sep <- exomePeakCalling(merip_bams = aln,
#'                         txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                         bsgenome = Hsapiens)
#'
#' sep <- normalizeGC(sep)
#'
#' sep <- glmDM(sep)
#'
#'
#' @name glmDM
#'
#' @rdname glmDM
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

  stopifnot( ( any(sep$design_Treatment) & any(!sep$design_Treatment) ) )

  if(glm_type == "auto") {
    if( all( table(colData(sep)$design_IP,colData(sep)$design_Treatment) > 1 ) ) {
      glm_type <- "DESeq2"
    } else {
      glm_type <- "Poisson"
    }
  }

  if(glm_type == "Poisson") {
    message("Differential modification analysis with Poisson GLM...")
  }

  if(glm_type == "NB") {
    message("Differential modification analysis with NB GLM...")
  }

  if(glm_type == "DESeq2") {
    message("Differential modification analysis with DESeq2...")
  }

  if(is.null(colData( sep )$sizeFactor)) {

    sep <- estimateSeqDepth(sep)

  }

  indx_mod <- grepl("mod", rownames( sep ) )

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

  if (!is.null(GCsizeFactors( sep ))) {

   if(LFC_shrinkage == "none") {

     DS_result <- suppressMessages( results(dds) )

   } else {

     DS_result  <- suppressMessages( lfcShrink(dds = dds, coef = 4, type = LFC_shrinkage) )

   }

    quantification_rst <- matrix( NA, nrow = nrow(SE_M[indx_mod,]), ncol = ncol(DS_result) )

    colnames(quantification_rst) <- colnames(DS_result)

    quantification_rst <- as.data.frame(quantification_rst)

    quantification_rst[(!gc_na_indx)[indx_mod],] <- as.data.frame( DS_result )

  } else {

    if(LFC_shrinkage == "none") {

    quantification_rst  <- as.data.frame( results(dds) )

    } else {

    suppressWarnings( quantification_rst <- as.data.frame( lfcShrink(dds=dds, coef=4, type = LFC_shrinkage) ) )

    }

  }

  rownames( quantification_rst ) = rownames( SE_M )[indx_mod]

  DESeq2Results( sep ) = as.data.frame( quantification_rst )

  return( sep )

})
