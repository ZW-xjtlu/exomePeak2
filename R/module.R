## A function to identify unmodified background using Gaussian mixture model
classifyBackground <- function(se, gmm_cut = 5){
  #require(mclust)
  IP_count <- assay(se[,se$IP_input == "IP"])
  input_count <- assay(se[,se$IP_input == "input"])

  input_sum <- rowSums(input_count)
  IP_sum <- rowSums(IP_count)
  indx_ratio <- (input_sum >= gmm_cut) & (IP_sum >= gmm_cut)
  logMratio <- log(IP_sum[indx_ratio] / input_sum[indx_ratio])
  rm(input_sum, IP_sum)

  gmm_fit <- Mclust(logMratio, G = 2)
  rm(logMratio)

  bg_class <- gmm_fit$parameters$mean
  bg_indx <- gmm_fit$classification == names(bg_class)[which.min(bg_class)]
  rm(gmm_fit, bg_class)
  rm(IP_count, input_count)

  rowData(se) <- DataFrame(bg = FALSE)
  rowData(se)$bg[which(indx_ratio)[bg_indx]] <- TRUE
  return(se)
}

## A function to estimate sequencing depth size factor from background
estimateColumnFactors <- function(se){
  stopifnot(!is.null(rowData(se)$bg))
  se$sf <- apply(assay(se)[rowData(se)$bg,], 2, function(x) median(x[x>0]))
  return(se)
}

## A function to fit GC content bias predictors using Poisson regression
fitBiasCurves <- function(se,
                          gc.knots = seq(from=.4, to=.6, length=3),
                          gc.bk = c(0,1)){
  #require(speedglm)
  #require(splines)
  model_lst <- vector("list", ncol(se))
  for( i in seq_len(ncol(se)) ){
    count_i <- assay(se)[,i]
    model_Matrix_i <- data.frame(x=rowData(se)$gc,y=count_i)
    model_lst[[i]] <- speedglm(y ~ ns(x, knots=gc.knots, Boundary.knots=gc.bk),
                               family = poisson(link="log"),
                               data = model_Matrix_i,
                               surface = "direct")
  }
  metadata(se) <- list(gc_model = model_lst,
                       fitm = metadata(se)[["fitm"]])
  return(se)
}

## A function to estimate feature specific size factors for GC content bias correction
estimateMatrixFactors <- function(se){
  #require(speedglm)
  #require(splines)
  if(!is.null(rowData(se)$gc)){
  sfm <- matrix(nrow = nrow(se), ncol = ncol(se))
  fitm <- matrix(nrow = 200, ncol = ncol(se))
  fit_range <- data.frame(x=seq(quantile(rowData(se)$gc, 0.05, na.rm = TRUE),
                                quantile(rowData(se)$gc, 0.95, na.rm = TRUE),
                                length.out = 200))

  for( i in seq_len(ncol(se)) ){
    fit_i <- metadata(se)[["gc_model"]][[i]]
    fitted_y <- predict.speedglm(fit_i, newdata = data.frame(x=rowData(se)$gc))
    sfm[,i] <- exp(fitted_y)/mean(exp(fitted_y))  # 0 center offsets
    fitted_p <- predict.speedglm(fit_i, newdata = fit_range) #fitted values to generate plot
    fitm[,i] <- exp(fitted_p)/mean(exp(fitted_p))
  }
  metadata(se) <- list(gc_model = metadata(se)[["gc_model"]],
                       fitm = fitm)
  }else{
    sfm <- matrix(1, nrow = nrow(se), ncol = ncol(se))
  }
  sfm <- t(t(sfm) * se$sf) #incorporate sequencing depth size factors
  assays(se, withDimnames=FALSE)$sfm <- sfm

  return(se)
}

## A function to count reads overlapped with features
featuresCounts <- function(features,
                           bam_dirs,
                           strandness = c("unstrand",
                                          "1st_strand",
                                          "2nd_strand"),
                           parallel = 1,
                           yield_size = 5000000){
  exist_indx <- file.exists(bam_dirs)

  if( any(!exist_indx) ){
    stop(paste0("Cannot find BAM file(s) of \n",
                paste(bam_dirs[!exist_indx],collapse = ", "),
                "\nplease re-check the input directories"))
  }

  #require(GenomicAlignments)
  #require(BiocParallel)
  strandness <- match.arg(strandness)

  ## Setup parallel number
  register(SerialParam())
  if (.Platform$OS.type == "windows") {
      register(SnowParam(workers = parallel))  # Windows alternative
  } else {
      suppressWarnings( register(MulticoreParam(workers = parallel)) )  # Unix/macOS
  }

  ## Setup bam file list
  bam_lst = BamFileList(file = bam_dirs,
                        asMates = TRUE)
  yieldSize(bam_lst) = yield_size

  ## Count using summarizeOverlaps (= HTSeq count Union)
  if(strandness == "1st_strand"){
    preprocess_func <- readsReverse
  }else{
    preprocess_func <- NULL
  }
  se <- summarizeOverlaps(
    features = features,
    reads = bam_lst,
    mode = "Union",
    inter.feature = FALSE,
    singleEnd = FALSE,
    preprocess.reads = preprocess_func,
    ignore.strand = strandness == "unstrand",
    fragments = TRUE
  )
  return(se)
}

#A function to calculate GC contents in fragment range
calculateGCcontents <- function(se,
                                fragment_length = 100,
                                exByGene,
                                genome){
  if(!is.null(genome)){
  rowData(se)$gc <- NA
  bin_size <- median(sum(width(rowRanges(se))))
  fragmentBins <- suppressWarnings(exonicFlank(rowRanges(se), exByGene, fragment_length))
  gcFreq <- letterFrequency(DNAStringSet(Views(genome, unlist(fragmentBins))), letters = "GC")
  gcFreq_gr <- tapply(gcFreq, rep(seq_along(fragmentBins), elementNROWS(fragmentBins)), sum)
  rowData(se)$gc[match(names(fragmentBins), rownames(se))] <- gcFreq_gr/sum(width(fragmentBins))
  se <- se[!is.na(rowData(se)$gc),]
  }
 return(se)
}

#A function to call peaks with GLM statistical testing
callPeaks <- function(se,
                      txdb,
                      test_method,
                      p_cutoff,
                      exByGene,
                      bin_size,
                      motif_based,
                      confounding_factor) {
  #Handling additional covariates
  if (!is.null(confounding_factor)) {
    if (is.factor(confounding_factor)) {
      confounding_factor <-
        data.frame(X = confounding_factor) #convert factor into single column data.frame
    }
    Dim_factor <- nrow(confounding_factor)
    if (Dim_factor != ncol(se)) {
      stop(
        "The provided `confounding_factor` variable has incompatible dimension with the total number of MeRIP-Seq samples (IP+input. Please double check.)"
      )
    }
    #set reference level with the design of confounding variables
    colData(se) <- cbind(DataFrame(confounding_factor), colData(se))
    se$IP_input <- relevel(factor(se$IP_input), "input")
    confounding_vars <- colnames(confounding_factor)
    formula_str <- paste("~", paste(confounding_vars, collapse = " + "), "+ IP_input")
    dds <- DESeqDataSet(se, as.formula(formula_str))
  }else{
    #set reference level
    se$IP_input <- relevel(factor(se$IP_input), "input")
    dds <- DESeqDataSet(se, ~ IP_input) 
  }

    #specify size factors
    if (is.null(assays(se)[["sfm"]])) {
      dds$sizeFactor <- se$sf
    } else{
      normalizationFactors(dds) <- assays(se)[["sfm"]]
    }

    #estimate dispersions
    if (test_method == "DESeq2") {
      dds <- estimateDispersions(dds)
    } else{
      dispersions(dds) <- 0
    }

    dds <- nbinomWaldTest(dds)
    res <- results(dds, altHypothesis = "greater")
    pvalue <- res$pvalue
    rm(dds)

    #Merge bins and annotate the resulting peaks
    pvalue[is.na(pvalue)] <- 1
    num_pos <- sum(pvalue < p_cutoff)
    peak <- reducePeaks(rowRanges(se)[pvalue <= p_cutoff], exByGene)
    if(!motif_based) peak <- peak[sum(width(peak)) >= bin_size]
    exbg <- exonsBy(txdb, by = "gene")

    mcols(peak) <- makePeakAnnot(peak, se, res, exbg)
    return(peak)
  }

#A function to call differential peaks with interactive GLM
callDiff <- function(se,
                     txdb,
                     test_method,
                     p_cutoff,
                     exByGene,
                     bin_size,
                     alt_hypothesis,
                     lfc_threshold,
                     motif_based,
                     absolute_diff,
                     confounding_factor){
  #Handling additional covariates
  if (!is.null(confounding_factor)) {
    if (is.factor(confounding_factor)) {
      confounding_factor <-
        data.frame(X = confounding_factor) #convert factor into single column data.frame
    }
    Dim_factor <- nrow(confounding_factor)
    if (Dim_factor != ncol(se)) {
      stop(
        "The provided `confounding_factor` variable has incompatible dimension with the total number of MeRIP-Seq samples (IP+input. Please double check.)"
      )
    }
    #set reference level with the design of confounding variables
    colData(se) <- cbind(DataFrame(confounding_factor), colData(se))
    se$IP_input <- relevel(factor(se$IP_input),"input")
    se$Perturbation <- relevel(factor(se$Perturbation),"C")
    confounding_vars <- colnames(confounding_factor)
    if(!absolute_diff){
      formula_str <- paste("~", paste(confounding_vars, collapse = " + "), "+ IP_input * Perturbation")
      dds <- DESeqDataSet(se, as.formula(formula_str))
      normalizationFactors(dds) <- assays(se)[["sfm"]]
    }else{
      formula_str <- paste("~", paste(confounding_vars, collapse = " + "), "+ Perturbation")
      dds <- DESeqDataSet(se[,se$IP_input!="input"], as.formula(formula_str))
      normalizationFactors(dds) <- assays(se[,se$IP_input!="input"])[["sfm"]]
    }
  }else{
    #Set reference levels
    se$IP_input <- relevel(factor(se$IP_input),"input")
    se$Perturbation <- relevel(factor(se$Perturbation),"C")
    if(!absolute_diff){
      dds <- DESeqDataSet(se, ~ IP_input * Perturbation)
      normalizationFactors(dds) <- assays(se)[["sfm"]]
    }else{
      dds <- DESeqDataSet(se[,se$IP_input!="input"], ~ Perturbation)
      normalizationFactors(dds) <- assays(se[,se$IP_input!="input"])[["sfm"]]
    }
  }
  
  #Fit differential models
  if(test_method == "DESeq2"){
    dds <- estimateDispersions(dds)
  }else{
    dispersions(dds) <- 0
  }

  dds <- nbinomWaldTest(dds)
  res <- results(dds,
                 altHypothesis = alt_hypothesis,
                 lfcThreshold = lfc_threshold)
  pvalue <- res$pvalue
  lfc <- res$log2FoldChange
  rm(dds)

  #Merge bins and annotate the resulting peaks
  pvalue[is.na(pvalue)] <- 1
  peak <- reducePeaks(rowRanges(se)[pvalue <= p_cutoff], exByGene)
  if(!motif_based) peak <- peak[sum(width(peak)) >= bin_size]
  exbg <- exonsBy(txdb, by = "gene")

  mcols(peak) <- makePeakAnnot(peak, se, res, exbg)
  return(peak)
}
