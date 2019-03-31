#' @title Exome peak calling on MeRIP-seq datasets while considering the biological variabilities.
#'
#' @description \code{exomePeakCalling} call RNA modification peaks on exome regions with statistical tests that account for the biological variabilities between samples.
#'
#' @details The function conduct exome level peak calling based on the read alignment results and the transcript annotations.
#'
#' @param merip_bams a \code{MeripBamFileList} object.
#' @param mode a character for the scope of peak calling on genome, can be one in "exon", "full_transcript", and "whole_genome".
#'
#' If mode == "whole_genome", the size of the genome is recommended to be less than 1e7.
#'
#' @param txdb a \code{TxDb} object, it can also be a single character string such as "hg19" which is processed by \code{\link{makeTxDbFromUCSC}}.
#' @param bsgenome a \code{\link{BSgenome}} object for the genome sequence, alternatively it could be the name of the reference genome recognized by \code{\link{getBSgenom}}.
#' @param gene_anno a string, which specifies a gene annotation GFF/GTF file if available, default: NA.
#' @param fragment_length a positive integer of the expected fragment length in bp; default 100.
#' @param binding_length a positive integer of the antibody binding length in IP samples; default 25.
#' @param step_length a positive integer of the shift size of the sliding window; default is the binding length.
#' @param glm_type a character, which can be one of the "DESeq2", "poisson", "NB". This argument specify the type of generalized linear model used in peak calling; Default to be "DESeq2".
#'
#' Under the default setting, the DESeq2 method is implemented on experiments with > 1 biological replicates for both IP and input samples.
#' The poisson GLM will be implemented otherwise.
#'
#' @param pc_count_cutoff a non negative integer value of the minimum average reads count per window used in peak calling; default 5.
#' @param gc_count_cutoff a non negative integer value of the minimum average reads count per window used in GC effect estimation; default 50.
#' @param p_cutoff a value of the p value cut-off used in peak calling; default NULL.
#' @param p_adj_cutoff a value of the adjusted p value cutoff used in DESeq inference; default 0.05.
#' @param logFC_cutoff a non negative numeric value of the log2 fold change (log2 IP/input) cutoff used in the inferene of peaks.
#' @param peak_width positive integer of the minimum width for the merged peaks; default \code{fragment_length} .
#' @param parallel a logical indicating whether to use parallel computation, consider this if your computer has more than 16GB RAM.
#' @param mod_annot a \code{GRanges} object for user provided single based RNA modification annotation. If provided, the peak calling step will be skipped.
#' Reads count will be performed using the provided annotation flanked by length of floor(fragment_length - binding_length/2).
#'
#' @param manual_background a \code{GRanges} object for user provided RNA modification background.
#' @param m6Aseq_background a logical of whether to use the topology knowledge of m6A-Seq to find appropriate background in GC effect estimation;
#' If TRUE, the GC effect will be estimated on bins that is not overlapping with long exons (exon length >= 400bp) and 5'UTR.
#' Also, the returned background control ranges returned will also exclude those regions;
#' It should not be select if you are analyzing MeRIP-seq data of other modification, such as hm5C; default TRUE.
#'
#' The background regions used in this senario will be the disjoint exon regions of the flanked provided sites.
#'
#' @param background a \code{GRanges} or \code{GRangesList} object for user provided background control regions on the genome.
#'
#' @return This function will return a \code{SummarizedExomePeak} object storing the ranges and reads counts information of the called peaks.
#'
#' @examples
#'
#'
#' @seealso \code{\link{exomePeakCalling}}
#' @import GenomicAlignments
#' @importFrom Rsamtools asMates
#' @import GenomicRanges
#' @import BiocParallel
#' @import SummarizedExperiment
#' @docType methods
#'
#' @name exomePeakCalling
#'
#' @rdname exomePeakCalling
#'
#' @export
#'
setMethod("exomePeakCalling",
          "MeripBamFileList",
          function(merip_bams = NULL,
                   txdb = NULL,
                   bsgenome = NULL,
                   glm_type = c("DESeq2", "NB", "poisson"),
                   background = c("Gaussian_mixture", "m6Aseq_prior", "manual", "all"),
                   manual_background = NULL,
                   gene_annot = NULL,
                   mod_annot = NULL,
                   fragment_length = 100,
                   binding_length = 25,
                   step_length = binding_length,
                   pc_count_cutoff = 5,
                   gc_count_cutoff = 50,
                   p_cutoff = NULL,
                   p_adj_cutoff = 0.05,
                   logFC_cutoff = 0,
                   peak_width = fragment_length / 2,
                   parallel = FALSE,
                   bp_param = NULL
          ) {

            ######################################################
            #                  Parameter check                   #
            ######################################################


            glm_type <- match.arg(glm_type)

            stopifnot(fragment_length > 0)

            stopifnot(step_length > 0)

            stopifnot(peak_width > 0)

            stopifnot(logFC_cutoff >= 0)

            stopifnot(pc_count_cutoff >= 0)

            stopifnot(gc_count_cutoff >= 0)

            if(!is.null(mod_annot)){
            stopifnot(is(mod_annot,"GRanges")|is(mod_annot,"GRangesList"))
            }

            if (is.null(bsgenome)) {
              warning(
                "Missing bsgenome argument, peak calling without GC content correction...",
                call. = FALSE,
                immediate. = TRUE
              )
            }

            if (!is.null(gene_annot)) {
              txdb <- makeTxDbFromGFF(gene_annot)
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
              colData(SE_Peak_counts) = DataFrame(metadata(merip_bams))


              #GC
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

              #Bin width
              rowData(SE_Peak_counts)$region_widths <-
                sum(width(rowRanges(SE_Peak_counts)))

              #Count index
              rowData(SE_Peak_counts)$indx_gc_est <-
                rowMeans(assay(SE_Peak_counts)) >= gc_count_cutoff


              ######################################################
              #                 Background search                  #
              ######################################################


              #Model based clustering

              m6A_prior = F

              if (background == "Gaussian_mixture") {
                message("Identifying background with Gaussian mixture model...")

                #suppress output...
                tc <- textConnection(NULL, "w")
                sink(tc)
                rowData(SE_Peak_counts)$indx_bg <- mclust_bg(se_peak_counts = SE_Peak_counts)
                sink()
                close(tc)

                if (sum(rowData(SE_Peak_counts)$indx_gc_est &
                        rowData(SE_Peak_counts)$indx_bg) < 30) {
                  warning(
                    "Background bin # < 30 using mclust, search background with m6A-seq prior",
                    call. = FALSE,
                    immediate. = FALSE
                  )
                   m6A_prior = T
                }

              }


              #Prior knowledge on m6A topology

              #Check for minimum background #


              if (background == "m6Aseq_prior" | m6A_prior) {
                indx_UTR5 <-
                  rowRanges(SE_Peak_counts) %over% fiveUTRsByTranscript(txdb)

                indx_longexon <-
                  rowRanges(SE_Peak_counts) %over% exons(txdb)[width(exons(txdb)) >= 400]

                rowData(SE_Peak_counts)$indx_bg = !(indx_UTR5 | indx_longexon)

                rm(indx_UTR5, indx_longexon, m6A_prior)
              }

              #Provided granges

              if (background == "manual") {

                rowData(SE_Peak_counts)$indx_bg = rowRanges(SE_Peak_counts) %over% manual_background

              }

              if (background == "all") {

                rowData(SE_Peak_counts)$indx_bg = T

              }

              #Check for minimum background #
              if (sum(rowData(SE_Peak_counts)$indx_gc_est &
                      rowData(SE_Peak_counts)$indx_bg) < 30) {
                warning(
                  "background bin # < 30 using m6A-seq prior, peak calling without background",
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

              if (glm_type == "poisson") {
                message("Peak calling with poisson GLM...")
              }

              if (glm_type == "NB") {
                message("Peak calling with NB GLM...")
              }

              if (glm_type == "DESeq2") {
                if (any(table(SE_Peak_counts$design_IP) == 1)) {
                  warning(
                    "At least one of the IP or input samples have no replicate, peak calling method changed to poisson GLM.",
                    call. = FALSE,
                    immediate. = TRUE
                  )
                  glm_type = "poisson"
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
                logFC_cutoff = logFC_cutoff,
                txdb = txdb
              )

              #Filter peak by width
              grl_mod <- grl_mod[sum(width(grl_mod)) >= peak_width]

              ######################################################
              #                 Round 2 reads count                #
              ######################################################

              #Flank the peaks

              gr_mod_flanked <- flank_on_exons(
                grl = grl_mod,
                flank_length = fragment_length - binding_length,
                txdb = txdb,
                index_flank = FALSE
              )

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
                }
              }

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

              mod_annot_flanked <- flank_on_exons(
                grl = mod_annot,
                flank_length = floor(fragment_length - binding_length /
                                       2),
                txdb = txdb,
                index_flank = FALSE
              )

              mod_annot_count <- disj_background(
                mod_gr = mod_annot_flanked,
                txdb = txdb,
                cut_off_num = 30,
                background_bins = rowRanges(SE_Peak_counts)[rowData(SE_Peak_counts)$indx_bg, ],
                background_types = background,
                control_width = peak_width
              )

              rm(SE_Peak_counts,mod_annot_flanked)

              message("Counting reads using single based annotation on exons...")

              if (!parallel) {
                register(SerialParam())
                register(MulticoreParam(workers = 1))
                register(SnowParam(workers = 1))
              } else {
                if (!is.null(bp_param)) {
                  register(bp_param, default = TRUE)
                }
              }

              SE_temp <- summarizeOverlaps(
                features = mod_annot_count,
                reads = merip_bams,
                param = Parameter(merip_bams),
                mode = "Union",
                inter.feature = FALSE,
                preprocess.reads = ifelse((LibraryType(merip_bams) == "1st_strand"), reads_five_POS_rev, reads_five_POS),
                singleEnd = !any(asMates(merip_bams)),
                ignore.strand = (LibraryType(merip_bams) == "unstranded"),
                fragments = any(asMates(merip_bams))
              )

              #Replace the rowRanges with the single based GRangesList.
              index_mod <- grepl("mod_", rownames(SE_temp))

              index_sb_annot <-
                as.numeric(gsub("mod_", "", rownames(SE_temp)[index_mod]))

              mod_count <-
                assay(SE_temp)[index_mod, ][match(as.numeric(names(mod_annot)), index_sb_annot), ]

              rownames(mod_count) <-
                paste0("mod_", seq_len(nrow(mod_count)))

              #Fill the rows that are not on exons with zero counts.(what it for?)
              mod_count[is.na(mod_count)] <- 0

              control_count <- assay(SE_temp)[!index_mod, ]

              #replace the metadata collumn with gene_ids
              mod_annot_gr <- unlist(mod_annot)

              mcols(mod_annot_gr) <- DataFrame(gene_id = NA)

              mcols(mod_annot_gr)$gene_id[index_sb_annot] <-
                mcols(unlist(rowRanges(SE_temp)))$gene_id[index_sb_annot]

              mod_annot <-
                split(mod_annot_gr, names(mod_annot_gr))

              mod_annot <-
                mod_annot[order(as.numeric(names(mod_annot)))]

              names(mod_annot) <- paste0("mod_", names(mod_annot))

              SummarizedExomePeaks <- SummarizedExperiment(
                assay = rbind(mod_count, control_count),

                rowRanges = c(mod_annot,
                              rowRanges(SE_temp)[!index_mod]),

                colData = DataFrame(metadata(merip_bams))

              )

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
