#' @title Peak Calling and Differential Analysis of MeRIP-seq.
#'
#' @description \code{exomePeak2} conducts peak calling and differential methylation analysis using BAM files of aligned MeRIP-seq reads.
#'
#' @details \code{\link{exomePeak2}} call (differential) RNA modification peaks and calculate peak statistics from \strong{BAM} files of a MeRIP-seq experiment.
#'
#' The transcript annotation (from either the \code{\link{TxDb}} object or the \strong{GFF} file) should be provided to perform analysis on exons.
#'
#' The genome name or \code{\link{BSgenome}} object is required to perform the GC content bias correction. If the \code{genome} argument is not provided (\code{= NULL}), the analysis will proceed without GC correction.
#'
#' If the \strong{BAM} files in treated samples are provided at the arguments \code{bam_ip_treated} and \code{bam_input_treated}, the statistics of differential modification detection on peaks/sites will be reported.
#'
#' Under the default setting, \code{\link{exomePeak2}} will save the results of (differential) modification analysis under a folder named \code{'exomePeak2_output'}.
#' The results generated include a \strong{BED} file, a \strong{RDS} file, and a \strong{CSV} table that stores the locations and statistics of the (differential) modified peaks/sites.
#'
#' @param bam_ip a \code{character} vector for the BAM file directories of IP samples.
#' @param bam_input a \code{character} vector for the BAM file directories of input samples.
#' @param bam_ip_treated a \code{character} vector for the BAM file directories of treated IP samples.
#' @param bam_input_treated a \code{character} vector for the BAM file directories of treated input samples.
#'
#' The arguments of \code{BAM_ip} and \code{BAM_input} are only required for differential methylation analysis.
#'
#' @param strandness a \code{character} specifying the strand protocol type of the RNA-seq library, can be one of \code{c("unstrand", "1st_strand", "2nd_strand")}; default \code{= "unstrand"}.
#'
#' \describe{
#' \item{\strong{unstrand}}{The randomly primed RNA-seq library type, i.e. both the strands generated during the first and the second strand sythesis are sequenced; example: Standard Illumina.}
#' \item{\strong{1st_strand}}{The first strand-specific RNA-seq library, only the strand generated during the first strand sythesis is sequenced; examples: dUTP, NSR, NNSR.}
#' \item{\strong{2nd_strand}}{The second strand-specific RNA-seq library, only the strand generated during the second strand sythesis is sequenced; examples: Ligation, Standard SOLiD.}
#' }
#'
#' @param txdb a \code{\link{TxDb}} object for the transcript annotation.
#'
#' @param genome a \code{character} or a \code{\link{BSgenome}} for the reference genome.
#'
#' The character should be the UCSC genome name which is acceptable by \code{\link{getBSgenome}} or/and \code{\link{makeTxDbFromUCSC}}; example: \code{"hg19"}.
#'
#' @param gff optional, a \code{character} which specifies the directory toward a gene annotation GFF/GTF file, it is applied when the \code{TxDb} object is not available; default \code{= NULL}.
#'
#' @param fragment_length a positive integer number for the expected fragment length (in bp); default \code{= 100}.
#'
#' @param bin_size a positive integer number for the width of the sliding window; default \code{= 25}.
#'
#' @param step_size a positive integer number for the step size of the sliding window; default \code{= 25}.
#'
#' @param test_method a \code{character} for the statistical testing method used in peak calling and differential analysis, can be one of c("Poisson", "DESeq2"); Default \code{= "Poisson"}
#'
#' \describe{
#' \item{\strong{\code{Poisson}}}{Wald test of Poisson GLM.}
#'
#' \item{\strong{\code{DESeq2}}}{Wald test of negative bionomial GLM with regularized estimation of over-dispersion parameters (implemented by DESeq2),}
#' }
#'
#' Note that when using the test method of DESeq2, a larger p-value cut-off (e.g. 0.001) is often required. The cutoff can be set via the argument \code{p_cutoff}.
#'
#' @param p_cutoff a \code{numeric} value for the p value cutoff in peak calling; default \code{= 1e-10}.
#'
#' @param parallel a \code{numeric} value specifying the number of cores for parallel computing; default \code{= 1}.
#'
#' @param plot_gc a \code{logical} for saving the plots of bins' GC content v.s. bins' fitted coverage curves, which can be used as a diagnosis for GC content bias; default \code{= TRUE}.
#'
#' @param save_output a \code{logical} for saving the outcomes on disk; default \code{= TRUE}.
#'
#' @param save_dir a \code{character} for the output directory; default \code{= getwd}.
#'
#' @param experiment_name a \code{character} for the folder name generated in the output directory that contains all the results; default: \code{="exomePeak2_output"}
#'
#' @param mode a \code{character} specifies the scope of peak calling on genome, can be one of \code{c("exon", "full_transcript", "whole_genome")}; Default \code{= "exon"}.
#'
#' \describe{
#' \item{\strong{\code{exon}}}{generate sliding windows over exonic regions.}
#'
#' \item{\strong{\code{full_transcript}}}{generate sliding windows over the full transcripts (include both introns and exons).}
#'
#' \item{\strong{\code{whole_genome}}}{generate sliding windows over the whole genome (include introns, exons, and the intergenic regions).}
#' }
#'
#' P.S. The full transcript mode and the whole genome mode demand big memory size (> 4GB) for large genomes.
#'
#' @param motif_based a \code{logical} for detecting (differential) modification over sites of motif; default \code{= FALSE}.
#' If \code{ = TRUE}, sliding windows will be replaced into the single based sites of the modification motif.
#'
#' @param motif_sequence a \code{character} for the motif sequence used for the reference sites, it is only applied when \code{motif_based = TRUE}; default \code{= "DRACH"}.
#' 
#' @param absolute_diff a \code{logical} for performing absolute differential modification without normalization over input control samples.
#' If \code{ = TRUE}, the regression design for differential modification test will be changed into comparing the direct changes of IP samples between treatment and control conditions; default \code{= FALSE}. 
#'
#' @return
#' a \code{\link{GRangesList}} object, the statistics and other annotations are saved in its metadata columns, which can be accessed through \code{mcol()}.
#' If \code{save_output = TRUE}, exomePeak2 will output results both as BED, CSV, and RDS files on disk.
#'
#' @examples
#'
#' ## Specify File Directories
#' GENE_ANNO_GTF = system.file("extdata", "example.gtf", package="exomePeak2")
#'
#' f1 = system.file("extdata", "IP1.bam", package="exomePeak2")
#' f2 = system.file("extdata", "IP2.bam", package="exomePeak2")
#' f3 = system.file("extdata", "IP3.bam", package="exomePeak2")
#' f4 = system.file("extdata", "IP4.bam", package="exomePeak2")
#' IP_BAM = c(f1,f2,f3,f4)
#' f1 = system.file("extdata", "Input1.bam", package="exomePeak2")
#' f2 = system.file("extdata", "Input2.bam", package="exomePeak2")
#' f3 = system.file("extdata", "Input3.bam", package="exomePeak2")
#' INPUT_BAM = c(f1,f2,f3)
#'
#' ## Peak Calling
#' res <- exomePeak2(bam_ip = IP_BAM,
#'                   bam_input = INPUT_BAM,
#'                   gff = GENE_ANNO_GTF,
#'                   genome = "hg19")
#' res        ## Peak ranges
#' mcols(res) ## Peak statistics
#'
#' ## Differential Peak Detection (Comparison of Two Conditions)
#'
#' f1 = system.file("extdata", "treated_IP1.bam", package="exomePeak2")
#' TREATED_IP_BAM = c(f1)
#' f1 = system.file("extdata", "treated_Input1.bam", package="exomePeak2")
#' TREATED_INPUT_BAM = c(f1)
#'
#' res <- exomePeak2(bam_ip = IP_BAM,
#'                   bam_input = INPUT_BAM,
#'                   bam_ip_treated = TREATED_IP_BAM,
#'                   bam_input_treated = TREATED_INPUT_BAM,
#'                   gff = GENE_ANNO_GTF,
#'                   genome = "hg19")
#' res        ## Peak ranges
#' mcols(res) ## Peak statistics
#'
#'
#' @name exomePeak2
#' @rdname exomePeak2
#'
#' @importFrom S4Vectors DataFrame na.omit queryHits subjectHits metadata runValue
#' @import GenomicFeatures
#' @import GenomicRanges
#' @importFrom BiocGenerics order
#' @import DESeq2
#' @importFrom magrittr "%>%"
#' @import BSgenome
#' @import BiocParallel
#' @import ggplot2
#' @import GenomicAlignments
#' @import mclust
#' @import speedglm
#' @import splines
#' @import SummarizedExperiment
#' @importFrom Biostrings DNAStringSet
#' @importFrom GenomeInfoDb seqlengths seqlevelsStyle "seqlengths<-" 
#' @importFrom IRanges subsetByOverlaps
#' @importFrom methods as is new
#' @importFrom stats relevel quantile median poisson
#' @importFrom utils capture.output read.table write.csv
#' @import Rsamtools
#' @importFrom rtracklayer export
#'
#' @export
#'
exomePeak2 <- function(bam_ip = NULL,
                       bam_input = NULL,
                       bam_ip_treated = NULL,
                       bam_input_treated = NULL,
                       txdb = NULL,
                       genome = NULL,
                       gff = NULL,
                       strandness = c("unstrand","1st_strand","2nd_strand"),
                       fragment_length = 100,
                       bin_size = 25,
                       step_size= 25,
                       test_method = c("Poisson", "DESeq2"),
                       p_cutoff = 1e-10,
                       parallel = 1,
                       plot_gc = TRUE,
                       save_output = TRUE,
                       save_dir = getwd(),
                       experiment_name = "exomePeak2_output",
                       mode = c("exon","full_transcript","whole_genome"),
                       motif_based = FALSE,
                       motif_sequence = "DRACH",
                       absolute_diff = FALSE
                      ){
  # Check input validity
  mode <- match.arg(mode)
  stopifnot(fragment_length > 0)
  stopifnot(step_size > 0)
  stopifnot(bin_size > 0)
  stopifnot(is.character(genome)|is(genome, "BSgenome")|is.null(genome))
  stopifnot(file.exists(save_dir))

  # Prepare transcript annotation
  if (is.null(gff) & is.null(txdb) & is.null(genome)){
    stop("Require one of the argument in txdb, gff, and genome for transcript annotation.")
  }

  if (is.null(txdb) & is.null(gff) & is.character(genome)) {
    txdb <- makeTxDbFromUCSC(genome)
  }

  if (is.null(txdb) & !is.null(gff)) {
    txdb <- makeTxDbFromGFF(gff) %>% quiet
  }

  if (mode != "exon") {
    txdb <- convertTxDb(txdb, type = mode) %>% quiet
  }

  # Prepare reference genome
  if (is.character(genome)){
    genome <- getBSgenome(genome)
  }

  # Run core function
  if (is.null(bam_ip_treated)) {
   res <- peakCalling(
      bam_IP = bam_ip,
      bam_input = bam_input,
      txdb = txdb,
      genome = genome,
      bin_size = bin_size,
      step_size = step_size,
      fragment_length = fragment_length,
      strandness = strandness,
      test_method = test_method,
      p_cutoff = p_cutoff,
      plot_gc = plot_gc,
      parallel = parallel,
      motif_based = motif_based,
      motif_sequence = "DRACH",
      fig_dir = file.path(save_dir, experiment_name),
      mode = mode
    )
   if(save_output) savePeak(res,
                            file.path(save_dir, experiment_name),
                            "peaks")
  } else {
    res <-  diffAnalysis(
        bam_IP = bam_ip,
        bam_input = bam_input,
        bam_IP_treated = bam_ip_treated,
        bam_input_treated = bam_input_treated,
        txdb = txdb,
        genome = genome,
        bin_size = bin_size,
        step_size = step_size,
        fragment_length = fragment_length,
        strandness = strandness,
        test_method = test_method,
        p_cutoff = p_cutoff,
        plot_gc = plot_gc,
        parallel = parallel,
        motif_based = motif_based,
        motif_sequence = "DRACH",
        absolute_diff = absolute_diff,
        fig_dir = file.path(save_dir, experiment_name),
        mode = mode
      )
    if(save_output) savePeak(res,
                             file.path(save_dir, experiment_name),
                             "diffPeaks")
  }

  return(res)
}
