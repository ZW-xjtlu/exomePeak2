## Differential peak calling with MeRIP-Seq.
## return the peaks in GRangesList format.
diffAnalysis <- function(bam_IP,
                         bam_input,
                         bam_IP_treated,
                         bam_input_treated,
                         txdb = NULL,
                         genome = NULL,
                         bin_size = 25,
                         step_size = 25,
                         fragment_length = 100,
                         strandness = c("unstrand", "1st_strand", "2nd_strand"),
                         gff = NULL,
                         test_method = c("Poisson", "DESeq2"),
                         p_cutoff = 1e-10,
                         lfc_threshold = 0,
                         diff_p_cutoff = 0.01,
                         alt_hypothesis = c("greaterAbs", "lessAbs", "greater", "less"),
                         plot_gc = FALSE,
                         parallel = 1,
                         motif_based = FALSE,
                         motif_sequence = "DRACH",
                         absolute_diff = FALSE,
                         fig_dir = "exomePeak2_output",
                         mode = c("exon","full_transcript","whole_genome")){
  #require(GenomicRanges)
  #require(GenomicFeatures)
  #require(SummarizedExperiment)
  #require(DESeq2)

  strandness <- match.arg(strandness)
  alt_hypothesis <- match.arg(alt_hypothesis)
  test_method <- match.arg(test_method)
  stopifnot(fragment_length > 0 & bin_size > 0 & step_size > 0)

  if(is.null(gff) & is.null(txdb)){
    stop("Please at least provide one among gff and TxDb for transcript annotation.")
  }

  if(is.null(txdb) & !is.null(gff)){
    txdb <- makeTxDbFromGFF(gff)
  }

  if(!(all(file.exists(bam_IP)) &
       all(file.exists(bam_input)) &
       all(file.exists(bam_IP_treated)) &
       all(file.exists(bam_input_treated)))){
    stop("At least one bam file directories provided cannot be found, please check it.")
  }

  if(is.null(genome) & motif_based){
    stop("Motif cannot be extracted without providing the 'genome' argument, please try to find the reference genome.")
  }

  if(is.null(genome)){
    warning("Reference genome not provided, GC content bias is left uncorrected.")
  }

  #Extract bins for count
  message("Extract bin features ... ", appendLF = F)
  exByGene  <- exonsByiGenes(txdb) %>% quiet
  if(!motif_based){
    peakBins <- exonicBins(exByGene, bin_size, step_size) %>% quiet
  }else{
    peakBins <- sampleSequence(motif_sequence, exByGene, genome) %>% quiet
    peakBins <-  resize(peakBins, 1, fix = "center") %>% quiet
  }
  mcols(peakBins) <- NULL
  if(is(peakBins, "GRanges")) peakBins <- split(peakBins, seq_along(peakBins)) %>% quiet
  message("OK")

  #Count the bam files
  message("Count reads on bin features ... ", appendLF = F)
  bam_dirs <- c(bam_IP, bam_IP_treated, bam_input, bam_input_treated)
  se <- featuresCounts(peakBins, bam_dirs, strandness, parallel) %>% quiet
  message("OK")

  #Annotate SummarizedExperiment
  length_indx <- c(length(bam_IP), length(bam_IP_treated),
                   length(bam_input), length(bam_input_treated))
  se$IP_input <- rep(c("IP", "IP", "input", "input"), length_indx)
  se$Perturbation <- rep(c("C", "Treated", "C", "Treated"), length_indx)
  rm(peakBins, length_indx)


  #Identify Backgrounds
  message("Identify background features ... ", appendLF = F)
  se <- classifyBackground(se) %>% quiet
  message("OK")

  #Estimate sample size factors
  message("Estimate sample sepecific size factors from the background ... ", appendLF = F)
  se <- estimateColumnFactors(se) %>% quiet
  message("OK")

  if(!is.null(genome)){
  #Calculate GC contents
  message("Calculate bin GC contents on exons ... ", appendLF = F)
  se <- calculateGCcontents(se, fragment_length, exByGene, genome) %>% quiet
  message("OK")

  #Fit GC content biases
  message("Fit GC curves with smoothing splines ... ", appendLF = F)
  se <- fitBiasCurves(se) %>% quiet
  message("OK")

  #Estimate matrix correction factors
  message("Calculate offset matrix for bins ... ", appendLF = F)
  se <- estimateMatrixFactors(se) %>% quiet
  message("OK")

  ## Plot GC bias fits
  if(plot_gc) plotGCbias(se, fig_dir) %>% quiet

  }else{
  #Assign matrix correction factors without GC offsets
  se <- estimateMatrixFactors(se) %>% quiet
  }
  
  #Filter low count rows if not exon mode
  if(mode %in% c("full_transcript","whole_genome")){
    se <- se[rowMeans(assay(se)) >= 5,]
  }
  
  #Peak calling
  message("Detect peaks with GLM ... ", appendLF = F)
  peaks <- callPeaks(se, txdb, test_method, p_cutoff, exByGene, bin_size, motif_based) %>% quiet
  message("OK")

  #Count the bam files
  message("Count reads on peaks ... ", appendLF = F)
  se2 <- featuresCounts(peaks, bam_dirs, strandness, parallel) %>% quiet
  rm(bam_dirs, peaks)
  message("OK")

  #Re-annotate SummarizedExperiment
  colData(se2) <- colData(se)
  metadata(se2) <- metadata(se)
  rm(se)

  #Calculate GC contents
  se2 <- calculateGCcontents(se2, fragment_length, exByGene, genome) %>% quiet

  #Estimate matrix correction factors
  message("Calculate offset matrix for peaks ... ", appendLF = F)
  se2 <- estimateMatrixFactors(se2) %>% quiet
  message("OK")

  #Differential calling
  if(!absolute_diff){
    message("Detect differentially modified peaks with interactive GLM ... ", appendLF = F)  
  }else{
    message("Detect absolute differentially modified peaks with GLM ... ", appendLF = F)  
  }
  diffPeaks <- callDiff(se2, txdb, test_method, diff_p_cutoff, exByGene, bin_size, alt_hypothesis, lfc_threshold, motif_based, absolute_diff) %>% quiet
  message("OK")

  return(diffPeaks)
}

