#' @title Export the methylation / differential methylation sites
#' @param sep A SummarizedExomePeak object.
#' @param format A character string of one of the c("tsv", "BED", "RDS"), indicating the exported format.
#'
#' - The choice tsv will save a tab separated values (tsv) file with the ranges and test information.
#'
#' - The choice BED will save a BEDGraph file with the score column being the -log2 adjusted p value.
#'
#' - The choice RDS will save a Rdata of the SummarizedExperiment object which will additional include a comprehensive summary of the count and the design information (recommended).
#'
#' @param save_dir The name of the file being saved; Default "exomePeak2_output".
#'
#' @param cut_off_pvalue A number between 0 and 1 indicate the p value cutoff in the exported result; Default NULL.
#'
#' @param cut_off_padj A number between 0 and 1 indicate the adjusted p value cutoff in the exported result; Default 0.05.
#'
#' @param cut_off_log2FC A non negative number indicating the log2 fold change cutoff of the exported result,
#'
#' only sites with log2 IP/input fold change bigger than this value are kept; Default 0.
#'
#' For differential methylation analysis, the absolute value of the log2 Odds ratio will be filtered by \code{cut_off_log2FC}.
#'
#' @param min_num_of_positive The minimum number of reported sites.
#' If the sites is filtered less than this number by its p values or effect sizes,
#' more sites will be reported by the order of the p value until it reaches this number.
#'
#' @param expected_direction The expected differential methylation direction, could be "hyper", "hypo", or "both".
#' This argument is useful when the treated group involves the perturbation of a writer or eraser protein for the modification; Default "both".
#'
#' @param inhibit_filter Remove all the filters on the result, this is desired for user provided modification annotation; Default FALSE.
#'
#' @param table_style Determine the style of the tsv table being exported, could be one of "bed" and "granges", the later would index the site containing multiple ranges with an id.
#'
#' @return a file containing the modification site genomic location, gene ids, statistics, and effect sizes will be created under the current working directory.
#'
#' @importFrom rtracklayer export
#' @import GenomicRanges
#'
#'@docType methods
#'
#'@name exportResults
#'
#'@rdname exportResults
#'
#' @export
#'

setMethod("exportResults",
          "SummarizedExomePeak",
          function(
                  sep,
                  format = c("tsv","BED","RDS"),
                  save_dir = "exomePeak2_output",
                  cut_off_pvalue = NULL,
                  cut_off_padj = 0.05,
                  cut_off_log2FC = 0,
                  min_num_of_positive = 2000,
                  expected_direction = "both",
                  inhibit_filter = FALSE,
                  table_style = c("bed","granges")
){

            if(!dir.exists(save_dir)) {
              dir.create(save_dir)
            }
            if (!any(sep$design_Treatment)) {
              file_name <- "meth"
            } else{
              file_name <- "diff_meth"
            }

            if (!inhibit_filter)
              file_name <- paste0("sig_", file_name)


            #In case of users have not run inference on the sep.
            if (is.null(DESeq2Results(sep))) {
              if (any(sep$design_Treatment)) {
                sep <- glm_dm(sep)
              } else{
                sep <- glm_meth(sep)
              }
            }

            #Check the validity of the arguments
            stopifnot(cut_off_pvalue <= 1 & cut_off_pvalue >= 0)
            stopifnot(cut_off_log2FC >= 0)
            format <- match.arg(format)
            table_style <- match.arg(table_style)

            #Decision for methylation

            if (!any(sep$design_Treatment)) {
              decision_meth <- decision_deseq2(
                Inf_RES = DESeq2Results(sep),
                log2FC_cut = cut_off_log2FC,
                P_cut = cut_off_pvalue,
                Padj_cut = cut_off_padj,
                Min_mod = min_num_of_positive
              )

              #In case of no sites are reported, export all the p values that are not NA

              index_keep <-
                which(
                  (DESeq2Results(sep)[[decision_meth$Cut_By_expected]] < decision_meth$Cut_Val_expected) &
                    (DESeq2Results(sep)$log2FoldChange > cut_off_log2FC)
                )

            } else {
              decision_dm <- decision_deseq2(
                Inf_RES = DESeq2Results(sep),
                log2FC_cut = cut_off_log2FC,
                P_cut = cut_off_pvalue,
                Padj_cut = cut_off_padj,
                Min_mod = min_num_of_positive,
                Exp_dir = expected_direction
              )

              if (expected_direction == "both") {
                indx_es <-
                  (abs(DESeq2Results(sep)$log2FoldChange) > cut_off_log2FC)
              } else {
                if (expected_direction == "hyper") {
                  indx_es <- (DESeq2Results(sep)$log2FoldChange > cut_off_log2FC)
                } else {
                  indx_es <- (DESeq2Results(sep)$log2FoldChange < -1 * cut_off_log2FC)
                }
              }

              index_keep <-
                which(DESeq2Results(sep)[[decision_dm$Cut_By_expected]] < decision_dm$Cut_Val_expected &
                        indx_es)

              if (length(index_keep) == 0) {
                stop(
                  "No sites could be left using the current filter, please change into a less rigorous one."
                )
              }

              if (inhibit_filter)
                index_keep <- rep(T, sum(grepl("meth_", rownames(sep))))

            }
            #now, create the final result summary that contain GRangesList with metadata collumns.
            result_grl <- rowRanges(sep)[grepl("meth_", rownames(sep))][index_keep]
            result_stat <- DESeq2Results(sep)[index_keep, ]


            #Make some clarifications on the final output results:
            #1. remove the .x in gene ids.
            #2. group names should be changed into mod names mod + id that is corresponding to the original table.
            #3. change log2FoldChange into log2 Odds ratio for differential methylation.
            result_gr <- unlist(result_grl)
            result_gr$gene_id <- gsub("\\.[0-9]*$", "", result_gr$gene_id)
            names(result_gr) <- gsub("\\..*$", "", names(result_gr))
            names(result_gr) <- gsub("meth", "mod", names(result_gr))

            if (!any(sep$design_Treatment)) {
              colnames(result_stat)[colnames(result_stat) == "log2FoldChange"] = "log2BetaValue"
              colnames(result_stat)[colnames(result_stat) == "lfcSE"] = "lbvSE"
              rownames(result_stat) = NULL

            } else {
              colnames(result_stat)[colnames(result_stat) == "log2FoldChange"] = "log2OddsRatio"
              colnames(result_stat)[colnames(result_stat) == "lfcSE"] = "lorSE"
              rownames(result_stat) = NULL

            }

            result_grl <- split(result_gr, names(result_gr))
            mcols(result_grl) <- result_stat

            if (format == "RDS") {
              result_se <- sep[grepl("meth_", rownames(sep))][index_keep]
              rowRanges(result_se) <- result_grl

              id_num <- as.numeric(gsub("^.*_", "", rownames(result_se)))
              id_index <- order(id_num)
              renamed_id <-
                paste0("mod_", rep(seq_along(id_num), table(id_num[id_index])))
              result_se <- result_se[id_index, ]
              rownames(result_se) <- renamed_id

              saveRDS(result_se, paste0(save_dir, "/", file_name, ".rds"))

            } else{
              id_num <- as.numeric(gsub("^.*_", "", names(result_grl)))
              id_index <- order(id_num)
              renamed_id <-
                paste0("mod_", rep(seq_along(id_num), table(id_num[id_index])))
              result_grl <- result_grl[id_index, ]
              names(result_grl) <- renamed_id

              #sort grl

              if (format == "tsv") {
                if (table_style == "granges") {
                  result_df <- as.data.frame(result_grl)
                  result_df <- result_df[, colnames(result_df) != "group"]
                  colnames(result_df)[colnames(result_df) == "group_name"] = "mod_name"
                  result_df <-
                    cbind(result_df, as.data.frame(mcols(result_grl))[rep(seq_along(result_grl), elementNROWS(result_grl)), ])

                } else {
                  scores <- -1 * log2(mcols(result_grl)$padj)

                  scores[is.na(scores)] <- 0

                  mcols(result_grl)$score <- scores

                  export(
                    object = result_grl,
                    con = paste0(save_dir, "/", file_name, ".bed"),
                    format = "BED"
                  )

                  result_df <-
                    read.table(paste0(save_dir, "/", file_name, ".bed"),
                               header = F,
                               sep = "\t")

                  colnames(result_df) <- c(
                    "chr",
                    "chromStart",
                    "chromEnd",
                    "name",
                    "score",
                    "strand",
                    "thickStart",
                    "thickEnd",
                    "itemRgb",
                    "blockCount",
                    "blockSizes",
                    "blockStarts"
                  )

                  mcols(result_grl) <-
                    mcols(result_grl)[, !colnames(mcols(result_grl)) %in% "score"]

                  result_df <- cbind(result_df , as.data.frame(mcols(result_grl)))

                }


                write.table(
                  result_df,
                  file = paste0(save_dir, "/", file_name, ".txt"),
                  sep = "\t",
                  row.names = FALSE,
                  col.names = TRUE
                )

              } else {
                scores <- -1 * log2(mcols(result_grl)$padj)

                scores[is.na(scores)] <- 0

                mcols(result_grl)$score <- scores

                export(
                  object = result_grl,
                  con = paste0(save_dir, ".", tolower(format)),
                  format = format
                )

              }
            }
          }
)
