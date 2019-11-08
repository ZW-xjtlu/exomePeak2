#' @title Perform Peak Calling on MeRIP-seq Dataset.
#'
#' @description \code{exomePeakCalling} call peaks of RNA modification from a MeRIP-seq data set.
#'
#' @details \code{exomePeakCalling} perform peak calling from the MeRIP-seq BAM files on exon regions defined by the user provided transcript annotations.
#' If the \code{\link{BSgenome}} object is provided, the peak calling will be conducted with the GC content bias correction.
#'
#' Under the default setting, for each window, exomePeak2 will fit a GLM of Negative Binomial (NB) with regulated estimation of the overdispersion parameters developed in \code{\link{DESeq}}.
#' Wald tests with H0 of IP/input Log2 Fold Change (LFC) <= 0 are performed on each of the sliding windows.
#' The significantly modified peaks are selected using the cutoff of p value < 0.0001.
#'
#' @param merip_bams a \code{MeripBamFileList} object returned by \code{\link{scanMeripBAM}}.
#'
#' @param txdb a \code{\link{TxDb}} object for the transcript annotation,
#' If the \code{TxDb} object is not available, it could be a \code{character} string of the UCSC genome name which is acceptable by \code{\link{makeTxDbFromUCSC}}, example: \code{"hg19"}.
#'
#' @param bsgenome a \code{\link{BSgenome}} object for the genome sequence information,
#' If the \code{BSgenome} object is not available, it could be a \code{character} string of the UCSC genome name which is acceptable by \code{\link{getBSgenome}}, example: \code{"hg19"}.
#'
#' @param gff_dir optional, a \code{character} which specifies the directory toward a gene annotation GFF/GTF file, it is applied when the \code{TxDb} object is not available, default \code{= NULL}.
#'
#' @param fragment_length a positive integer number for the expected fragment length in nucleotides; default \code{= 100}.
#'
#' @param binding_length a positive integer number for the expected binding length of the anti-modification antibody in IP samples; default \code{= 25}.
#'
#' @param step_length a positive integer number for the shift distances of the sliding window; default \code{= binding_length}.
#'
#' @param glm_type a \code{character} speciefies the type of Generalized Linear Model (GLM) fitted for the purpose of statistical inference during peak calling, which can be one of the \code{c("DESeq2", "NB", "Poisson")}.
#'
#' \describe{
#' \item{\strong{\code{DESeq2}}}{Fit the GLM defined in function \code{\link{DESeq}}, which is the NB GLM with regulated estimation of the overdispersion parameters.}
#'
#' \item{\strong{\code{NB}}}{Fit the Negative Binomial (NB) GLM.}
#'
#' \item{\strong{\code{Poisson}}}{Fit the Poisson GLM.}
#' }
#'
#' By default, the DESeq2 GLMs are fitted on the data set with > 1 biological replicates for both the IP and input samples, the Poisson GLM will be fitted otherwise.
#'
#' @param pc_count_cutoff a \code{numeric} value for the cutoff on average window's reads count in peak calling; default \code{= 5}.
#'
#' @param bg_count_cutoff a \code{numeric} value for the cutoff on average window's reads count in background identification; default \code{= 50}.
#'
#' @param p_cutoff a \code{numeric} value for the cutoff on p values in peak calling; default \code{= 0.0001}.
#'
#' @param p_adj_cutoff a \code{numeric} value for the cutoff on Benjamini Hochberg adjusted p values in peak calling; default \code{= NULL}.
#'
#' @param log2FC_cutoff a \code{numeric} value for the cutoff on log2 IP over input fold changes in peak calling; default \code{= 1}.
#'
#' @param consistent_peak a \code{logical} of whether the positive peaks returned should be consistent among all the replicates; default \code{= FALSE}.
#'
#' @param consistent_log2FC_cutoff a \code{numeric} for the modification log2 fold changes cutoff in the peak consisency calculation; default = 1.
#'
#' @param consistent_fdr_cutoff a \code{numeric} for the BH adjusted C-test p values cutoff in the peak consistency calculation; default { = 0.05}. Check \code{\link{ctest}}.
#'
#' @param alpha a \code{numeric} for the binomial quantile used in the consitent peak filter; default\code{ = 0.05}.
#'
#' @param p0 a \code{numeric} for the binomial proportion parameter used in the consistent peak filter; default \code{= 0.8}.
#'
#' For a peak to be consistently methylated, the minimum number of significant enriched replicate pairs is defined as the 1 - alpha quantile of a binomial distribution with p = p0 and N = number of possible pairs between replicates.
#'
#' The consistency defined in this way is equivalent to the rejection of an exact binomial test with null hypothesis of p < p0 and N = replicates number of IP * replicates number of input.
#'
#' @param peak_width a \code{numeric} value for the minimum width of the merged peaks; default \code{= fragment_length} .
#'
#' @param parallel a \code{logical} value of whether to use parallel computation, typlically it requires more than 16GB of RAM if \code{parallel = TRUE}; default \code{= FALSE}.
#'
#' @param mod_annot a \code{\link{GRanges}} object for user provided single based RNA modification annotation.
#'
#' If user provides the single based RNA modification annotation, this function will perform reads count on the provided annotation flanked by length \code{= floor(fragment_length - binding_length/2)}.
#'
#' @param manual_background  a \code{\link{GRanges}} object for the user provided unmodified background; default \code{= NULL}.
#'
#' @param correct_GC_bg a \code{logical} value of whether to estimate the GC content linear effect on background regions; default \code{= FALSE}.
#'
#' If \code{= TRUE}, it could lead to a more accurate estimation of GC content bias for the RNA modifications that are highly biologically related to GC content.
#'
#' @param qtnorm a \code{logical} of whether to perform subset quantile normalization after the GC content linear effect correction; default \code{= FASLE}.
#'
#' If \code{qtnorm = TRUE}, subset quantile normalization will be applied within the IP and input samples seperately to account for the inherent differences between the marginal distributions of IP and input samples.
#'
#' @param background a \code{character} specifies the method for the background finding, i.e. to identify the windows without modification signal. It could be one of \code{c("Gaussian_mixture", "m6Aseq_prior", "manual", "all")};  default \code{= "Gaussian_mixture"}.
#'
#' In order to accurately account for the technical variations, it is often neccessary to estimate the GC content linear effects on windows without modification signals (background).
#'
#' The following methods are supported in \code{ExomePeak2} to differentiate the no modification background windows from the modification containig windows.
#'
#' \describe{
#'  \item{\strong{\code{Gaussian_mixture}}}{The background is identified by Multivariate Gaussian Mixture Model (MGMM) with 2 mixing components on the vectors containing methylation level estimates and GC content, the background regions are predicted by the Bayes Classifier on the learned GMM.}
#'
#'  \item{\strong{\code{m6Aseq_prior}}}{The background is identified by the prior knowledge of m6A topology, the windows that are not overlapped with long exons (exon length >= 400bp) and 5'UTR are treated as the background windows.
#'
#'  This type of background should not be used if the MeRIP-seq data is not using anti-m6A antibody.
#'
#'  }
#'
#'  \item{\strong{\code{manual}}}{The background regions are defined by user manually at the argument \code{manual_background}.}
#'
#'  \item{\strong{\code{all}}}{Use all windows as the background. This is equivalent to not differentiating background and signal.
#'  It can lead to biases during the sequencing depth and the GC content correction factors estimation.
#'  }
#' }
#'
#' @param bp_param optional, a \code{\link{BiocParallelParam}} object that stores the configuration parameters for the parallel execution.
#'
#' @return a \code{\link{SummarizedExomePeak}} object.
#'
#' @examples
#'
#' # Load packages for the genome sequence and transcript annotation
#'
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#'
#' # Peak Calling
#'
#' exomePeakCalling(
#'   merip_bams = scanMeripBAM(
#'   bam_ip = c("IP_rep1.bam",
#'              "IP_rep2.bam",
#'              "IP_rep3.bam"),
#'   bam_input = c("input_rep1.bam",
#'                 "input_rep2.bam",
#'                 "input_rep3.bam"),
#'   paired_end = TRUE),
#'   txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'   bsgenome = Hsapiens)
#'
#' # Analysis with single based modification annotation
#'
#' annot_dir <- system.file("extdata", "m6A_hg19_annot.rds", package = "exomePeak2")
#'
#' m6A_hg19_gr <- readRDS(annot_dir)
#'
#' exomePeakCalling(
#'   merip_bams = scanMeripBAM(
#'   bam_ip = c("IP_rep1.bam",
#'              "IP_rep2.bam",
#'              "IP_rep3.bam"),
#'   bam_input = c("input_rep1.bam",
#'                 "input_rep2.bam",
#'                 "input_rep3.bam"),
#'  paired_end = TRUE),
#'  txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'  bsgenome = Hsapiens,
#'  mod_annot = m6A_hg19_gr)
#'
#'
#' @seealso \code{\link{exomePeak2}}, \code{\link{glmM}}, \code{\link{glmDM}}, \code{\link{normalizeGC}}, \code{\link{exportResults}}
#' @import GenomicAlignments
#' @importFrom Rsamtools asMates
#' @import GenomicRanges
#' @import BiocParallel
#' @import SummarizedExperiment
#'
#' @aliases exomePeakCalling
#'
#' @rdname exomePeakCalling-methods
#'
#' @export
#'
setMethod("exomePeakCalling",
          "MeripBamFileList",
          function(merip_bams = NULL,
                   txdb = NULL,
                   bsgenome = NULL,
                   mod_annot = NULL,
                   glm_type = c("DESeq2", "NB", "Poisson"),
                   background = c("all",
                                  "Gaussian_mixture",
                                  "m6Aseq_prior",
                                  "manual"),
                   manual_background = NULL,
                   correct_GC_bg = FALSE,
                   qtnorm = FALSE,
                   gff_dir = NULL,
                   fragment_length = 100,
                   binding_length = 25,
                   step_length = binding_length,
                   peak_width = fragment_length / 2,
                   pc_count_cutoff = 5,
                   bg_count_cutoff = 50,
                   p_cutoff = 0.0001,
                   p_adj_cutoff = NULL,
                   log2FC_cutoff = 1,
                   consistent_peak = FALSE,
                   consistent_log2FC_cutoff = 1,
                   consistent_fdr_cutoff = 0.05,
                   alpha = 0.05,
                   p0 = 0.8,
                   parallel = FALSE,
                   bp_param = NULL
          ) {

            ######################################################
            #                  Parameter check                   #
            ######################################################


            glm_type <- match.arg(glm_type)

            background <- match.arg(background)

            stopifnot(fragment_length > 0)

            stopifnot(step_length > 0)

            stopifnot(peak_width > 0)

            stopifnot(log2FC_cutoff >= 0)

            stopifnot(pc_count_cutoff >= 0)

            stopifnot(bg_count_cutoff >= 0)

            if(!is.null(mod_annot)){
            stopifnot(is(mod_annot,"GRanges")|is(mod_annot,"GRangesList"))
            }

            if (is.null(bsgenome)) {
              warning(
                "Missing BSgenome, peak calling without GC content correction...",
                call. = FALSE,
                immediate. = TRUE
              )
            }

            if (!is.null(gff_dir)) {
              txdb <- makeTxDbFromGFF(gff_dir)
            } else {
              if (is.null(txdb)) {
                stop("Missing transcript annotation, please provide TxDb object or GFF/GTF file...")
              }

              if (!is(txdb, "TxDb")) {
                txdb <- makeTxDbFromUCSC(txdb)
              }
            }

              message("Generating bins on exons...")


              ######################################################
              #             Sliding window generation              #
              ######################################################


              exome_bins_grl <- exome_bins_from_txdb(
                txdb = txdb,
                window_size = binding_length,
                step_size = step_length
              )

              message("Counting reads on bins...")


              ######################################################
              #                Initial reads count                 #
              ######################################################

              if (!parallel) {
                register(SerialParam())
                register(MulticoreParam(workers = 1))
                register(SnowParam(workers = 1))
              } else {
                if (!is.null(bp_param)) {
                  register(bp_param, default = TRUE)
                }
              }
              yieldSize(merip_bams) = 5000000 #Control the yield size to inhibit memory overflow

              SE_Peak_counts <- suppressWarnings( summarizeOverlaps(
                features = split_by_name(
                  flank_on_exons(
                    grl = exome_bins_grl,
                    flank_length = fragment_length - binding_length,
                    txdb = txdb,
                    index_flank = FALSE
                  )
                ),
                reads = merip_bams,
                param = Parameter(merip_bams),
                mode = "Union",
                inter.feature = FALSE,
                preprocess.reads = ifelse((LibraryType(merip_bams) == "1st_strand"), reads_five_POS_rev, reads_five_POS),
                singleEnd = !any(asMates(merip_bams)),
                ignore.strand = (LibraryType(merip_bams) == "unstranded"),
                fragments = any(asMates(merip_bams))
              ) )

              ######################################################
              #               Reads count annotation               #
              ######################################################

              #Experimental design
              colData(SE_Peak_counts) <- DataFrame(metadata(merip_bams))


              #GC content calculation
              if (is.null(bsgenome)) {

              } else{
                flanked_gr <- unlist(rowRanges(SE_Peak_counts))
                names(flanked_gr) <-
                  gsub("\\..*$", "", names(flanked_gr))
                GC_freq <-
                  as.vector(letterFrequency(Views(bsgenome, flanked_gr), letters = "CG"))
                sum_freq <- tapply(GC_freq, names(flanked_gr), sum)
                sum_freq <- sum_freq[names(rowRanges(SE_Peak_counts))]
                rowData(SE_Peak_counts)$gc_contents <-
                  sum_freq / sum(width(rowRanges(SE_Peak_counts)))
                rm(flanked_gr, GC_freq, sum_freq)
              }

              #Bin width calculation
              rowData(SE_Peak_counts)$region_widths <-
                sum(width(rowRanges(SE_Peak_counts)))

              #Count cutoff index
              rowData(SE_Peak_counts)$indx_gc_est <-
                rowMeans(assay(SE_Peak_counts)) >= bg_count_cutoff


              ######################################################
              #                 Background search                  #
              ######################################################

              #Model based clustering
              m6A_prior = F

              if (background == "Gaussian_mixture") {
                message("Identifying background with Gaussian Mixture Model...")

                rowData(SE_Peak_counts)$indx_bg <- quiet( mclust_bg(se_peak_counts = SE_Peak_counts) )

                if (sum(rowData(SE_Peak_counts)$indx_gc_est &
                        rowData(SE_Peak_counts)$indx_bg) < 30) {

                  warning(
                    "Number of the background bins < 30. Background method changed to the uniform.",
                    call. = FALSE,
                    immediate. = FALSE
                  )
                  rowData(SE_Peak_counts)$indx_bg = T
                }

              }


              #Define background using prior knowledge of m6A topology

              if (background == "m6Aseq_prior" | m6A_prior) {
                indx_UTR5 <-
                  rowRanges(SE_Peak_counts) %over% fiveUTRsByTranscript(txdb)

                indx_longexon <-
                  rowRanges(SE_Peak_counts) %over% exons(txdb)[width(exons(txdb)) >= 400]

                rowData(SE_Peak_counts)$indx_bg = !(indx_UTR5 | indx_longexon)

                rm(indx_UTR5, indx_longexon, m6A_prior)
              }

              #Define background using user provided GRanges

              if (background == "manual") {

                rowData(SE_Peak_counts)$indx_bg = rowRanges(SE_Peak_counts) %over% manual_background

              }

              if (background == "all") {

                rowData(SE_Peak_counts)$indx_bg = T

              }

              #Check for minimum background number
              if (sum(rowData(SE_Peak_counts)$indx_gc_est &
                      rowData(SE_Peak_counts)$indx_bg) < 30) {
                warning(
                  "Number of the background bins < 30. Background method changed to the uniform.\n",
                  call. = FALSE,
                  immediate. = FALSE
                )
                rowData(SE_Peak_counts)$indx_bg = T
              }


              if (is.null(mod_annot)) {

              ######################################################
              #                    Peak calling                    #
              ######################################################

              #Change bins into initial widths
              rowData_tmp <- rowData(SE_Peak_counts)

              rowRanges(SE_Peak_counts) <-
                exome_bins_grl[as.numeric(rownames(SE_Peak_counts))]

              rowData(SE_Peak_counts) <- rowData_tmp

              rm(exome_bins_grl,rowData_tmp)

              if (glm_type == "Poisson") {
                message("Peak calling with Poisson GLM...")
              }

              if (glm_type == "NB") {
                message("Peak calling with NB GLM...")
              }

              if (glm_type == "DESeq2") {
                if (any(table(SE_Peak_counts$design_IP) == 1)) {
                  warning(
                    "At least one of the IP or input group has no replicates. Peak calling method changed to Poisson GLM.",
                    call. = FALSE,
                    immediate. = TRUE
                  )
                  glm_type = "Poisson"
                } else{
                  message("Peak calling with DESeq2...")
                }
              }

              grl_mod <- call_peaks_with_GLM(
                SE_bins = SE_Peak_counts,
                glm_type = glm_type,
                count_cutoff = pc_count_cutoff,
                p_cutoff = p_cutoff,
                p_adj_cutoff = p_adj_cutoff,
                log2FC_cutoff = log2FC_cutoff,
                correct_GC_bg = correct_GC_bg,
                qtnorm =  qtnorm,
                txdb = txdb,
                consistent_peak = consistent_peak,
                consistent_log2FC_cutoff = consistent_log2FC_cutoff,
                consistent_fdr_cutoff = consistent_fdr_cutoff,
                alpha = alpha,
                p0 = p0
              )

              #Filter peak by width
              grl_mod <- grl_mod[sum(width(grl_mod)) >= peak_width]

              ######################################################
              #                 Round 2 reads count                #
              ######################################################

              #Flank the peaks

              gr_mod_flanked <- suppressWarnings( flank_on_exons(
                grl = grl_mod,
                flank_length = fragment_length - binding_length,
                txdb = txdb,
                index_flank = FALSE
              ) )

              #Set control regions with background disjoint by peaks

              count_row_features <- disj_background(
                mod_gr = gr_mod_flanked,
                txdb = txdb,
                cut_off_num = 30,
                background_bins = rowRanges(SE_Peak_counts)[rowData(SE_Peak_counts)$indx_bg, ],
                background_types = background,
                control_width = peak_width
              )

              rm(SE_Peak_counts, gr_mod_flanked)

              message("Counting reads on the merged peaks and the control regions...")

              if (!parallel) {
                register(SerialParam(), default = FALSE)
                register(MulticoreParam(workers = 1))
                register(SnowParam(workers = 1))
              } else {
                if (!is.null(bp_param)) {
                  register(bp_param, default = TRUE)
                } else {
                  register(SerialParam(), default = FALSE)
                  register(MulticoreParam(workers = 3))
                  register(SnowParam(workers = 3))
                }
              }

              yieldSize(merip_bams) = 5000000 #Control the yield size to inhibit memory overflow

              SummarizedExomePeaks <- suppressWarnings( summarizeOverlaps(
                features = count_row_features,
                reads = merip_bams,
                param = Parameter(merip_bams),
                mode = "Union",
                inter.feature = FALSE,
                preprocess.reads = ifelse((LibraryType(merip_bams) == "1st_strand"), reads_five_POS_rev, reads_five_POS),
                singleEnd = !any(asMates(merip_bams)),
                ignore.strand = (LibraryType(merip_bams) == "unstranded"),
                fragments = any(asMates(merip_bams))
              ) )

              #retrieve the set of unflanked modification sites to replace the row ranges.

              names(grl_mod) <- paste0("mod_", names(grl_mod))

              rowRanges(SummarizedExomePeaks)[seq_along(grl_mod)] <- grl_mod

              rm(count_row_features, grl_mod)

              colData(SummarizedExomePeaks) <- DataFrame(metadata(merip_bams))

              colnames(SummarizedExomePeaks) <- names(merip_bams)

            } else {

              ######################################################
              #             Count reads on annotation              #
              ######################################################

              stopifnot(is(mod_annot, "GRanges") |
                          is(mod_annot, "GRangesList"))

              stopifnot(all(width(mod_annot)) == 1)


              if (is(mod_annot, "GRanges")) {
                mod_annot <-
                  split(mod_annot , seq_along(mod_annot))

              } else {
                names(mod_annot) <-
                  seq_along(mod_annot) #make sure the annotation is indexed by integer sequence

              }

              mod_annot_flanked <- suppressWarnings( flank_on_exons(
                grl = mod_annot,
                flank_length = floor(fragment_length - binding_length / 2),
                txdb = txdb,
                index_flank = FALSE
              ) )

              mod_annot_flanked <- split_by_name(mod_annot_flanked)

              names(mod_annot_flanked) <- paste0("mod_", names(mod_annot_flanked))

              # mod_annot_count <- suppressWarnings( disj_background(
              #   mod_gr = mod_annot_flanked,
              #   txdb = txdb,
              #   cut_off_num = 30,
              #   background_bins = rowRanges(SE_Peak_counts)[rowData(SE_Peak_counts)$indx_bg, ],
              #   background_types = background,
              #   control_width = peak_width
              # ) )

              SE_Peak_counts_bg <- SE_Peak_counts[rowData(SE_Peak_counts)$indx_bg, ]

              rm(SE_Peak_counts)

              SE_Peak_counts_bg <- SE_Peak_counts_bg[rowMeans(assay(SE_Peak_counts_bg)[,colData(SE_Peak_counts_bg)$design_IP]) >= 50 &
                                                       rowMeans(assay(SE_Peak_counts_bg)[,!colData(SE_Peak_counts_bg)$design_IP]) >= 50 ]

              SE_Peak_counts_bg <- SE_Peak_counts_bg[!rowRanges(SE_Peak_counts_bg)%over%mod_annot, ]

              rownames(SE_Peak_counts_bg) <- paste0("control_", seq_len(dim(SE_Peak_counts_bg)[1]))

              message("Counting reads on the single based annotation...")

              if (!parallel) {
                register(SerialParam())
                register(MulticoreParam(workers = 1))
                register(SnowParam(workers = 1))
              } else {
                if (!is.null(bp_param)) {
                  register(bp_param, default = TRUE)
                } else {
                  register(SerialParam(), default = FALSE)
                  register(MulticoreParam(workers = 3))
                  register(SnowParam(workers = 3))
                }
              }

              yieldSize(merip_bams) = 5000000 #Control the yield size to inhibit memory overflow

              SE_temp <- summarizeOverlaps(
                features = c(mod_annot_flanked,rowRanges(SE_Peak_counts_bg)),
                reads = merip_bams,
                param = Parameter(merip_bams),
                mode = "Union",
                inter.feature = FALSE,
                preprocess.reads = ifelse((LibraryType(merip_bams) == "1st_strand"), reads_five_POS_rev, reads_five_POS),
                singleEnd = !any(asMates(merip_bams)),
                ignore.strand = (LibraryType(merip_bams) == "unstranded"),
                fragments = any(asMates(merip_bams))
              )

              rm(SE_Peak_counts_bg)

              #Replace the rowRanges with the single based GRangesList.
              index_mod <- grepl("mod_", rownames(SE_temp))

              index_se <- match(as.numeric(names(mod_annot)),
                                as.numeric(gsub("mod_", "", rownames(SE_temp)[index_mod])))

              mod_count <-
                assay(SE_temp)[index_mod, ][index_se, ]

              rm(index_se)

              rownames(mod_count) <-
                paste0("mod_", seq_len(nrow(mod_count)))

              #Fill the rows that are not on exons with zero counts.
              mod_count[is.na(mod_count)] <- 0

              control_count <- assay(SE_temp)[!index_mod, ]

              #replace the metadata collumn with gene_ids
              mod_annot_gr <- unlist(mod_annot)

              mcols(mod_annot_gr) <- DataFrame(gene_id = NA)

              SE_temp_gr <- unlist(rowRanges(SE_temp)[index_mod])

              mcols(mod_annot_gr)$gene_id <-
                mcols(SE_temp_gr)$gene_id[match(as.numeric(names(mod_annot_gr)),
                                            as.numeric(gsub("mod_.*\\.", "", names(SE_temp_gr))))]

              rm(SE_temp_gr)

              mod_annot <- split(mod_annot_gr, names(mod_annot_gr))

              mod_annot <- mod_annot[order(as.numeric(names(mod_annot)))]

              names(mod_annot) <- paste0("mod_", names(mod_annot))

              SummarizedExomePeaks <- SummarizedExperiment(
                assay = rbind(mod_count, control_count),

                rowRanges = c(mod_annot,
                              rowRanges(SE_temp)[!index_mod]),

                colData = DataFrame(metadata(merip_bams))

              )

              colnames(SummarizedExomePeaks) <- names(merip_bams)

            }

            #Annotate GC content if BSgenome is provided
            if (!is.null(bsgenome)) {
              elementMetadata(SummarizedExomePeaks) <- GC_content_over_grl(
                bsgenome = bsgenome,
                txdb = txdb,
                grl = rowRanges(SummarizedExomePeaks),
                fragment_length = fragment_length,
                binding_length = binding_length,
                effective_GC = FALSE
              )
            }

            return(
              new(
                "SummarizedExomePeak",
                rowRanges = SummarizedExomePeaks@rowRanges,
                colData = SummarizedExomePeaks@colData,
                assays = SummarizedExomePeaks@assays,
                NAMES = SummarizedExomePeaks@NAMES,
                elementMetadata = SummarizedExomePeaks@elementMetadata,
                metadata = SummarizedExomePeaks@metadata,
                DESeq2Results = data.frame()
              )
            )
        }
)
