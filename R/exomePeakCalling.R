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
#' @param txdb a \code{\link{TxDb}} object for the transcript annotation.
#'
#' If the \code{TxDb} object is not available, it could be a \code{character} string of the UCSC genome name which is acceptable by \code{\link{makeTxDbFromUCSC}}. For example: \code{"hg19"}.
#'
#' @param bsgenome a \code{\link{BSgenome}} object for the genome sequence information.
#'
#' If the \code{BSgenome} object is not available, it could be a \code{character} string of the UCSC genome name which is acceptable by \code{\link{getBSgenome}}. For example: \code{"hg19"}.
#'
#' @param genome a \code{character} string of the UCSC genome name which is acceptable by \code{\link{getBSgenome}} or/and \code{\link{makeTxDbFromUCSC}}. For example: \code{"hg19"}.
#'
#' By default, the argument = NA, it should be provided when the \code{BSgenome} or/and the \code{TxDb} object are not available.
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
#' @param p_cutoff a \code{numeric} value for the cutoff on p values in peak calling; default \code{= 1e-05}.
#'
#' @param p_adj_cutoff a \code{numeric} value for the cutoff on Benjamini Hochberg adjusted p values in peak calling; default \code{= NULL}.
#'
#' @param log2FC_cutoff a \code{numeric} value for the cutoff on log2 IP over input fold changes in peak calling; default \code{= 0}.
#'
#' @param min_peak_width a \code{numeric} value for the minimum width of the merged peaks; default \code{= fragment_length/2} .
#'
#' @param min_peak_width a \code{numeric} value for the maximum width of the merged peaks; default \code{= fragment_length*10} .
#' 
#' @param parallel a \code{numeric} value specifying the number of cores used for parallel computing; default \code{= 3}.
#'
#' @param mod_annot a \code{\link{GRanges}} or \code{\link{GRangesList}} object for user provided single based RNA modification annotation, the widths of the ranged object should be all equal to 1.
#'
#' If user provides the single based RNA modification annotation, exomePeak2 will perform reads count and (differential) modification quantification on the provided annotation.
#'
#' The single base annotation will be flanked by length = floor(fragment_length - binding_length/2) to account for the fragment length of the sequencing library.
#'
#' @param manual_background  a \code{\link{GRanges}} object for the user provided unmodified background; default \code{= NULL}.
#'
#' @param correct_GC_bg a \code{logical} value of whether to estimate the GC content linear effect on background regions; default \code{= TRUE}.
#'
#' If \code{= TRUE}, it could lead to a more accurate estimation of GC content bias for the RNA modifications that are highly biologically related to GC content.
#'
#' @param qtnorm a \code{logical} of whether to perform subset quantile normalization after the GC content linear effect correction; default \code{= FASLE}.
#'
#' If \code{qtnorm = TRUE}, subset quantile normalization will be applied within the IP and input samples seperately to account for the inherent differences between the marginal distributions of IP and input samples.
#'
#' @param background_method a \code{character} specifies the method for the background finding, i.e. to identify the windows without modification signal. It could be one of the "Gaussian_mixture", "m6Aseq_prior", "manual", and "all";  default \code{= "all"}.
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
#' ### Define File Directories
#'
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
#' f1 = system.file("extdata", "treated_IP1.bam", package="exomePeak2")
#' TREATED_IP_BAM = c(f1)
#' f1 = system.file("extdata", "treated_Input1.bam", package="exomePeak2")
#' TREATED_INPUT_BAM = c(f1)
#'
#' ### Peak Calling
#'
#' MeRIP_Seq_Alignment <- scanMeripBAM(
#'                          bam_ip = IP_BAM,
#'                          bam_input = INPUT_BAM,
#'                          paired_end = FALSE
#'                         )
#'
#' sep <- exomePeakCalling(
#'             merip_bams = MeRIP_Seq_Alignment,
#'             gff_dir = GENE_ANNO_GTF,
#'             genome = "hg19"
#'             )
#'
#' sep <- normalizeGC(sep)
#'
#' sep <- glmM(sep)
#'
#' exportResults(sep)
#'
#' ### Differential Modification Analysis (Comparison of Two Conditions)
#'
#' MeRIP_Seq_Alignment <- scanMeripBAM(
#'                          bam_ip = IP_BAM,
#'                          bam_input = INPUT_BAM,
#'                          bam_treated_ip = TREATED_IP_BAM,
#'                          bam_treated_input = TREATED_INPUT_BAM,
#'                          paired_end = FALSE
#'                         )
#'
#' sep <- exomePeakCalling(
#'             merip_bams = MeRIP_Seq_Alignment,
#'             gff_dir = GENE_ANNO_GTF,
#'             genome = "hg19"
#'             )
#'
#' sep <- normalizeGC(sep)
#'
#' sep <- glmDM(sep)
#'
#' exportResults(sep)
#'
#' @seealso \code{\link{exomePeak2}}, \code{\link{glmM}}, \code{\link{glmDM}}, \code{\link{normalizeGC}}, \code{\link{exportResults}}
#' @import GenomicAlignments
#' @importFrom Rsamtools asMates yieldSize<-
#' @import GenomicRanges
#' @importFrom BiocParallel register SerialParam MulticoreParam SnowParam
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
                   genome = NA,
                   mod_annot = NULL,
                   glm_type = c("DESeq2",
                                "NB",
                                "Poisson"),
                   background_method = c("all",
                                         "Gaussian_mixture",
                                         "m6Aseq_prior",
                                         "manual"),
                   manual_background = NULL,
                   correct_GC_bg = TRUE,
                   qtnorm = FALSE,
                   gff_dir = NULL,
                   fragment_length = 100,
                   binding_length = 25,
                   step_length = binding_length,
                   min_peak_width = fragment_length/2,
                   max_peak_width = fragment_length*10,
                   pc_count_cutoff = 5,
                   bg_count_cutoff = 50,
                   p_cutoff = 1e-05,
                   p_adj_cutoff = NULL,
                   log2FC_cutoff = 0,
                   parallel = 1,
                   bp_param = NULL
          ) {

            ######################################################
            #                  Parameter check                   #
            ######################################################


            glm_type <- match.arg(glm_type)

            background_method <- match.arg(background_method)

            stopifnot(fragment_length > 0)

            stopifnot(step_length > 0)

            stopifnot(min_peak_width > 0)
            
            stopifnot(max_peak_width > 0)

            stopifnot(log2FC_cutoff >= 0)

            stopifnot(pc_count_cutoff >= 0)

            stopifnot(bg_count_cutoff >= 0)

            if(!is.null(mod_annot)){
            stopifnot(is(mod_annot,"GRanges")|is(mod_annot,"GRangesList"))
            }

            if(!is.na(genome)) {
              if(!is(bsgenome,"BSgenome")) bsgenome = genome
              if(!is(txdb,"TxDb") & is.null(gff_dir)) txdb = genome
            }

            if (!is.null(gff_dir) & is.null(txdb)) {
              message("Make the TxDb object ... ", appendLF = FALSE)
              txdb <- suppressMessages( makeTxDbFromGFF(gff_dir) )
              message("OK")
            } else {
              if (is.null(txdb)) {
                stop("Missing transcript annotation, please provide the genome name or the transcript annotation package/files.")
              }

              if (!is(txdb, "TxDb")) {
                txdb <- makeTxDbFromUCSC(txdb)
              }
            }

            if(is.character(bsgenome)) {
              bsgenome <- getBSgenome(bsgenome)
            }

            if (is.null(bsgenome)) {
              warning(
                "Missing BSgenome or UCSC genome name, peak calling without GC content correction.",
                call. = FALSE,
                immediate. = TRUE
              )
            }
              message("Generate bins on exons ... ", appendLF = F)


              ######################################################
              #             Sliding window generation              #
              ######################################################


              exome_bins_grl <- exome_bins_from_txdb(
                txdb = txdb,
                window_size = binding_length,
                step_size = step_length
              )

              message("OK")

              message("Count reads on bins ... ", appendLF = F)


              ######################################################
              #                Initial reads count                 #
              ######################################################

              if (!is.null(bp_param)) {
                register(bp_param, default = TRUE)
              } else {
                register(SerialParam())
                suppressWarnings( register(MulticoreParam(workers = parallel)) )
                register(SnowParam(workers = parallel))
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

              message("OK")

              ######################################################
              #               Reads count annotation               #
              ######################################################

              #Add Experimental design
              colData(SE_Peak_counts) <- DataFrame(metadata(merip_bams))


              #GC content calculation

              if (is.null(bsgenome)) {

              } else {
                message("Calculate GC contents on exons ... ", appendLF = F)

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

                message("OK")

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

              if (background_method == "Gaussian_mixture") {
                message("Identify background with Gaussian Mixture Model ... ", appendLF = F)

                rowData(SE_Peak_counts)$indx_bg <- quiet( mclust_bg(se_peak_counts = SE_Peak_counts) )

                message("OK")

                if (sum(rowData(SE_Peak_counts)$indx_gc_est &
                        rowData(SE_Peak_counts)$indx_bg) < 30) {

                  warning(
                    "Number of the background bins < 30. Background method changed to 'All'.",
                    call. = FALSE,
                    immediate. = FALSE
                  )
                  rowData(SE_Peak_counts)$indx_bg = rowMeans(assay(SE_Peak_counts)) >= bg_count_cutoff
                }

              }


              #Define background using prior knowledge of m6A topology

              if (background_method == "m6Aseq_prior" | m6A_prior) {
                indx_UTR5 <-
                  rowRanges(SE_Peak_counts) %over% fiveUTRsByTranscript(txdb)

                indx_longexon <-
                  rowRanges(SE_Peak_counts) %over% exons(txdb)[width(exons(txdb)) >= 400]

                rowData(SE_Peak_counts)$indx_bg = !(indx_UTR5 | indx_longexon)

                rm(indx_UTR5, indx_longexon, m6A_prior)
              }

              #Define background using user provided GRanges

              if (background_method == "manual") {

                rowData(SE_Peak_counts)$indx_bg = rowRanges(SE_Peak_counts) %over% manual_background

              }

              if (background_method == "all") {

                rowData(SE_Peak_counts)$indx_bg = rowMeans(assay(SE_Peak_counts)) >= bg_count_cutoff

              } else {

              #Check for minimum background number
              if (sum(rowData(SE_Peak_counts)$indx_gc_est &
                      rowData(SE_Peak_counts)$indx_bg) < 30) {
                warning(
                  'Number of the background bins < 30. Background method changed to "All".\n',
                  call. = FALSE,
                  immediate. = FALSE
                )
                rowData(SE_Peak_counts)$indx_bg = rowMeans(assay(SE_Peak_counts)) >= bg_count_cutoff
               }
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

              if (glm_type == "DESeq2") {
                if (any(table(SE_Peak_counts$design_IP) == 1)) {
                  warning(
                    "At least one of the IP or input group has no replicates. Peak calling method changed to Poisson GLM.\n",
                    call. = TRUE,
                    immediate. = FALSE
                  )
                  glm_type = "Poisson"
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
                txdb = txdb
              )

              if(length(grl_mod) == 0) stop("No modification peaks are detected. Please try smaller values of `p_cutoff`, e.x. 0.01.")

              #Filter peak by width
              indx_keep <- (sum(width(grl_mod)) >= min_peak_width) & (sum(width(grl_mod)) <= max_peak_width)
              grl_mod <- grl_mod[indx_keep]
              rm(indx_keep)
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
                background_bins = rowRanges(SE_Peak_counts)[rowData(SE_Peak_counts)$indx_bg, ],
                background_types = background_method,
                control_width = min_peak_width
              )
              
              rm(SE_Peak_counts, gr_mod_flanked)
              message("Count reads on the merged peaks and the control regions ... ", appendLF = F)

              if (!is.null(bp_param)) {
                register(bp_param, default = TRUE)
              } else {
                register(SerialParam())
                suppressWarnings( register(MulticoreParam(workers = parallel)) )
                register(SnowParam(workers = parallel))
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

              message("OK")

              #retrieve the set of unflanked modification sites to replace the row ranges.

              names(grl_mod) <- paste0("peak_", names(grl_mod))

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

                names(mod_annot) <- NULL

                mod_annot <-
                  split(mod_annot , seq_along(mod_annot))

              }

                names(mod_annot) <-
                  paste0("peak_",seq_along(mod_annot)) #make sure the annotation is indexed by integer sequence


              mod_annot_flanked <- suppressWarnings( flank_on_exons(
                grl = mod_annot,
                flank_length = floor(fragment_length - binding_length / 2),
                txdb = txdb,
                index_flank = FALSE
              ) )

              mod_annot_flanked <- split_by_name(mod_annot_flanked)

              # mod_annot_count <- suppressWarnings( disj_background(
              #   mod_gr = mod_annot_flanked,
              #   txdb = txdb,
              #   background_bins = rowRanges(SE_Peak_counts)[rowData(SE_Peak_counts)$indx_bg, ],
              #   background_types = background,
              #   control_width = peak_width
              # ) )

              SE_Peak_counts_bg <- SE_Peak_counts[rowData(SE_Peak_counts)$indx_bg, ]

              rm(SE_Peak_counts)

              SE_Peak_counts_bg <- SE_Peak_counts_bg[!rowRanges(SE_Peak_counts_bg)%over%mod_annot, ]

              rownames(SE_Peak_counts_bg) <- paste0("control_", seq_len(dim(SE_Peak_counts_bg)[1]))

              message("Count reads on the single based annotation ... ", appendLF = F)

              if (!is.null(bp_param)) {
                register(bp_param, default = TRUE)
              } else {
                register(SerialParam())
                suppressWarnings( register(MulticoreParam(workers = parallel)) )
                register(SnowParam(workers = parallel))
              }

              yieldSize(merip_bams) = 5000000 #Control the yield size to inhibit memory overflow

              SE_temp <- summarizeOverlaps(
                features = mod_annot_flanked,
                reads = merip_bams,
                param = Parameter(merip_bams),
                mode = "Union",
                inter.feature = FALSE,
                preprocess.reads = ifelse((LibraryType(merip_bams) == "1st_strand"), reads_five_POS_rev, reads_five_POS),
                singleEnd = !any(asMates(merip_bams)),
                ignore.strand = (LibraryType(merip_bams) == "unstranded"),
                fragments = any(asMates(merip_bams))
              )

              #rm(SE_Peak_counts_bg)
              rm(mod_annot_flanked)

              message("OK")

              #Replace the rowRanges with the single based GRangesList.

              mod_annot <- unlist(mod_annot)

              mod_annot$gene_id <- NA

              indx_modc <- as.numeric(gsub("peak_","",rownames(SE_temp)))

              mod_annot$gene_id[rep(indx_modc,elementNROWS(rowRanges(SE_temp)))] <- unlist(rowRanges( SE_temp ))$gene_id

              names(mod_annot) <- NULL

              mod_annot <- split(mod_annot, seq_along(mod_annot))

              names(mod_annot) <- paste0("peak_", names(mod_annot))

              SummarizedExomePeaks <- SummarizedExperiment(
                assays = matrix(0,nrow = nrow(SE_Peak_counts_bg)+length(mod_annot),
                                 ncol = ncol(SE_Peak_counts_bg)),

                rowRanges = c(mod_annot,rowRanges(SE_Peak_counts_bg)),

                colData = DataFrame(metadata(merip_bams))

              )

              indx_ctl <- grepl("control",rownames(SummarizedExomePeaks))

              assay(SummarizedExomePeaks,withDimnames=FALSE)[indx_ctl,] <- assay(SE_Peak_counts_bg)

              rm(SE_Peak_counts_bg, indx_ctl)

              assay(SummarizedExomePeaks,withDimnames=FALSE)[indx_modc,] <- assay(SE_temp)

              rm(SE_temp, indx_modc)

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
              SummarizedExomePeak(
                rowRanges = rowRanges(SummarizedExomePeaks),
                colData = colData(SummarizedExomePeaks),
                assays = Assays(assays(SummarizedExomePeaks)),
                NAMES = NULL,
                elementMetadata = DataFrame(matrix(vector(), nrow(SummarizedExomePeaks), 0)),
                metadata = metadata(SummarizedExomePeaks),
                exomePeak2Results = data.frame()
              )
            )
        }
)
