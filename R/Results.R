#' @title Export the modification / differential modification sites
#' @param sep A SummarizedExomePeak object.
#' @param cut_off_pvalue A number between 0 and 1 indicate the p value cutoff in the exported result; Default NULL.
#' @param cut_off_padj A number between 0 and 1 indicate the adjusted p value cutoff in the exported result; Default 0.05.
#' @param cut_off_log2FC A non negative number indicating the log2 fold change cutoff of the exported result,
#' only sites with log2 IP/input fold change bigger than this value are kept; Default 0.
#'
#' For differential modification analysis, the absolute value of the log2 Odds ratio will be filtered by \code{cut_off_log2FC}.
#'
#' @param min_num_of_positive The minimum number of reported sites.
#' If the sites is filtered less than this number by its p values or effect sizes,
#' more sites will be reported by the order of the p value until it reaches this number.
#'
#' @param expected_direction The expected differential modification direction, could be "hyper", "hypo", or "both".
#' This argument is useful when the treated group involves the perturbation of a writer or eraser protein for the modification; default "both".
#'
#' @param inhibit_filter Remove all the filters upon on the result; Default TRUE.
#'
#' @param table_style Determine the style of the tsv table being returned, could be one of "bed" and "granges", the later would index the site containing multiple ranges with an id.
#'
#' @return a data.frame containing the modification site genomic location, gene ids, statistics, and effect sizes will be returned.
#'
#' @importFrom rtracklayer export
#' @import GenomicRanges
#'
#' @name Results
#'
#' @rdname Results
#'
#' @export
#'

setMethod("Results",
          "SummarizedExomePeak",
          function(sep,
                   cut_off_pvalue = NULL,
                   cut_off_padj = 0.05,
                   cut_off_log2FC = 0,
                   min_num_of_positive = 30,
                   expected_direction = "both",
                   inhibit_filter = TRUE,
                   table_style = c("bed", "granges")
                   ) {

            table_style <- match.arg(table_style)

            if (table_style == "bed") {
              if (file.access(".") != 0)
                stop(
                  "The current working directory is not accessible by R,
                  please consider change the table_style into 'granges'."
                )
            }

            #In case of users have not run inference on the sep.
            if (is.null(DESeq2Results(sep))) {
              if (any(sep$design_Treatment)) {
                sep <- glmDM(sep)
              } else{
                sep <- glm_M(sep)
              }
            }

            #Check the validity of the arguments
            stopifnot(cut_off_pvalue <= 1 & cut_off_pvalue >= 0)
            stopifnot(cut_off_log2FC >= 0)
            table_style <- match.arg(table_style)

            #Decision for modification

            if (inhibit_filter){
              index_keep <- rep(T, sum(grepl("mod_", rownames(sep))))
            } else {

            if (!any(sep$design_Treatment)) {
              decision_mod <- decision_deseq2(
                Inf_RES = DESeq2Results(sep),
                log2FC_cut = cut_off_log2FC,
                P_cut = cut_off_pvalue,
                Padj_cut = cut_off_padj,
                Min_mod = min_num_of_positive
              )

              #In case of no sites are reported, export all the p values that are not NA

              index_keep <-
                which(
                  (DESeq2Results(sep)[[decision_mod$Cut_By_expected]] < decision_mod$Cut_Val_expected) &
                    (DESeq2Results(sep)$log2FoldChange > cut_off_log2FC)
                )

            } else {
              decision_dm <- decision_deseq2(
                Inf_RES = DESeq2Results(sep),
                log2FC_cut = cut_off_log2FC,
                P_cut = cut_off_pvalue,
                Padj_cut = cut_off_padj,
                Min_mod = min(min_num_of_positive,nrow(DESeq2Results(sep))),
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
            }
            }

            #now, create the final result summary that contain GRangesList with metadata collumns.
            result_grl <-
              rowRanges(sep)[grepl("mod_", rownames(sep))][index_keep]
            result_stat <- DESeq2Results(sep)[index_keep, ]


            #Make some clarifications on the final output results:
            #1. remove the .x in gene ids.
            #2. group names should be changed into mod names mod + id that is corresponding to the original table.
            #3. change log2FoldChange into log2 Odds ratio for differential modification.
            result_gr <- unlist(result_grl)
            result_gr$gene_id <-
              gsub("\\.[0-9]*$", "", result_gr$gene_id)
            names(result_gr) <-
              gsub("\\..*$", "", names(result_gr))

            if (!any(sep$design_Treatment)) {
              colnames(result_stat)[colnames(result_stat) == "log2FoldChange"] = "mod.log2.fc"
              rownames(result_stat) = NULL

            } else {
              colnames(result_stat)[colnames(result_stat) == "log2FoldChange"] = "diff.log2.fc"
              rownames(result_stat) = NULL

            }

            result_grl <- split(result_gr, names(result_gr))
            mcols(result_grl) <- result_stat


            id_num <-
              as.numeric(gsub("^.*_", "", names(result_grl)))
            id_index <- order(id_num)
            renamed_id <-
              paste0("mod_", rep(seq_along(id_num), table(id_num[id_index])))
            result_grl <- result_grl[id_index, ]
            names(result_grl) <- renamed_id

            #sort grl

            if (table_style == "granges") {
              result_df <- as.data.frame(result_grl)
              result_df <-
                result_df[, colnames(result_df) != "group"]
              colnames(result_df)[colnames(result_df) == "group_name"] = "mod_name"
              result_df <-
                cbind(result_df, as.data.frame(mcols(result_grl))[rep(seq_along(result_grl), elementNROWS(result_grl)), ])

            } else {
              scores <- -1 * log2(mcols(result_grl)$padj)

              scores[is.na(scores)] <- 0

              mcols(result_grl)$score <- scores

              export(object = result_grl,
                     con = "temp___1.bed",
                     format = "BED")

              result_df <-
                read.table("temp___1.bed", header = F, sep = "\t")

              file.remove("temp___1.bed")

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

              result_df$geneID <- sapply( result_grl , function(x) x$gene_id[1] )

              mcols(result_grl) <-
                mcols(result_grl)[, !colnames(mcols(result_grl)) %in% "score"]

              result_df <-
                cbind(result_df , as.data.frame(mcols(result_grl)))

            }

            #Attach gene ids at the end of the table

            return(result_df)

            })
