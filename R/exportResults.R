#' @title Export the (Differential) Modification Peaks/Sites and their associated LFC Statistics
#' @param sep a \code{\link{SummarizedExomePeak}} object.
#' @param format a \code{character} for the exported format, can be one in \code{c("tsv", "BED", "RDS")}.
#'
#' \describe{
#'  \item{\strong{\code{tsv}}}{
#'  export a tab separated values (tsv) table with the genomic location and LFC statistics.
#'  }
#'
#'  \item{\strong{\code{BED}}}{
#'  export a BEDGraph file with the score column = -log2(adjusted p value).
#'  }
#'
#'  \item{\strong{\code{RDS}}}{
#'  export the RDS file of the \code{\link{SummarizedExperiment}} object.
#'  }
#' }
#'
#' @param save_dir a \code{character} for the name of the directory being saved; Default \code{= "exomePeak2_output"}.
#'
#' @param cut_off_pvalue a \code{numeric} value for the p value cutoff in the exported result; Default \code{= NULL}.
#'
#' @param cut_off_padj a \code{numeric} value for the adjusted p value cutoff in the exported result; Default \code{= 0.05}.
#'
#' @param cut_off_log2FC a \code{numeric} value for the log2 fold change cutoff of the exported result,
#' only the sites with abs(LFC) larger than this value are kept; Default \code{= 0}.
#'
#' @param min_num_of_positive a \code{numeric} value for the minimum number of reported sites.
#' If the number of remaining sites is less than this number after filtering, additional sites will be reported by the increasing order of the p value to meet this number.
#'
#' @param expected_direction a \code{character} for the expected differential modification direction, could be one in \code{c("hyper", "hypo", "both")}.
#'
#' \describe{
#'  \item{\strong{\code{hyper}}}{
#'  only report the peaks/sites with interactive LFC > 0.
#'  }
#'
#'  \item{\strong{\code{hypo}}}{
#'  only report the peaks/sites with interactive LFC < 0.
#'  }
#'
#'  \item{\strong{\code{both}}}{
#'  report the peaks/sites in both directions.
#'  }
#' }
#'
#' This argument is useful when the treated group involves the perturbation of a known writer or eraser protein; Default "both".
#'
#' @param inhibit_filter a \code{logical} for whether to remove all the filters, this option is useful when quantification on single based site annotation; Default \code{= FALSE}.
#'
#' @param table_style a \code{character} for the style of the tsv table being exported, could be one in \code{c("bed","granges")}.
#'
#' \describe{
#'  \item{\strong{\code{bed}}}{
#'  the genomic locations in the table are represented by BEDgraph style.
#'  }
#'
#'  \item{\strong{\code{granges}}}{
#'  the genomic locations in the table are represented by GRanges style.
#'  }
#' }
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
#' paired_end = TRUE
#' )
#'
#' sep <- exomePeakCalling(merip_bams = aln,
#'                         txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                         bsgenome = Hsapiens)
#'
#' sep <- normalizeGC(sep)
#'
#' sep <- glmM(sep)
#'
#' exportResults(sep)
#'
#' @importFrom rtracklayer export
#' @import GenomicRanges
#'
#' @docType methods
#'
#' @name exportResults
#'
#' @rdname exportResults
#'
#' @export
#'

setMethod("exportResults",
          "SummarizedExomePeak",
          function(sep,
                   format = c("tsv", "BED", "RDS"),
                   save_dir = "exomePeak2_output",
                   cut_off_pvalue = NULL,
                   cut_off_padj = 0.05,
                   cut_off_log2FC = 0,
                   min_num_of_positive = 30,
                   expected_direction = "both",
                   inhibit_filter = FALSE,
                   table_style = c("bed", "granges")) {

            if (!dir.exists(save_dir)) {
              dir.create(save_dir)
            }

            if (!any(sep$design_Treatment)) {
              file_name <- "mod"
            } else{
              file_name <- "diff_mod"
            }

            if (!inhibit_filter)
              file_name <- paste0("sig_", file_name)


            #In case of users have not run inference on the sep.
            if (is.null(DESeq2Results(sep))) {
              if (any(sep$design_Treatment)) {
                sep <- glmDM(sep)
              } else{
                sep <- glmM(sep)
              }
            }

            #Check the validity of the arguments
            stopifnot(cut_off_pvalue <= 1 & cut_off_pvalue >= 0)
            stopifnot(cut_off_log2FC >= 0)
            format <- match.arg(format)
            table_style <- match.arg(table_style)

            #Decision for modification

            if (inhibit_filter) {
              index_keep <- rep(T, sum(grepl("mod_", rownames(sep))))
            } else {

            if (!any(sep$design_Treatment)) {
              decision_mod <- decision_deseq2(
                Inf_RES = DESeq2Results(sep),
                log2FC_cut = cut_off_log2FC,
                P_cut = cut_off_pvalue,
                Padj_cut = cut_off_padj,
                Min_mod = min(min_num_of_positive, nrow(DESeq2Results(sep)))
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
            }
            }
            #now, create the final result summary that contain GRangesList with metadata collumns.
            result_grl <-
              rowRanges(sep)[grepl("mod_", rownames(sep))][index_keep]
            result_stat <- DESeq2Results(sep)[index_keep,]


            #Make some clarifications on the final output results:
            #1. remove the .x in gene ids.
            #2. group names should be changed into mod names mod + id that is corresponding to the original table.
            #3. change log2FoldChange into log2 Odds ratio for differential modification.
            result_gr <- unlist(result_grl)
            result_gr$gene_id <-
              gsub("\\.[0-9]*$", "", result_gr$gene_id)
            names(result_gr) <- gsub("\\..*$", "", names(result_gr))

            if (!any(sep$design_Treatment)) {
              colnames(result_stat)[colnames(result_stat) == "log2FoldChange"] = "mod.log2.fc"
              rownames(result_stat) = NULL

            } else {
              colnames(result_stat)[colnames(result_stat) == "log2FoldChange"] = "diff.log2.fc"
              rownames(result_stat) = NULL

            }

            result_grl <- split(result_gr, names(result_gr))
            mcols(result_grl) <- result_stat

            if (format == "RDS") {
              result_se <- sep[grepl("mod_", rownames(sep))][index_keep]
              rowRanges(result_se) <- result_grl

              id_num <-
                as.numeric(gsub("^.*_", "", rownames(result_se)))
              id_index <- order(id_num)
              renamed_id <-
                paste0("mod_", rep(seq_along(id_num), table(id_num[id_index])))
              result_se <- result_se[id_index,]
              rownames(result_se) <- renamed_id


              saveRDS(result_se, file.path(save_dir, paste0(file_name, ".rds")))


            } else{
              id_num <- as.numeric(gsub("^.*_", "", names(result_grl)))
              id_index <- order(id_num)
              renamed_id <-
                paste0("mod_", rep(seq_along(id_num), table(id_num[id_index])))
              result_grl <- result_grl[id_index,]
              names(result_grl) <- renamed_id

              #sort grl

              if (format == "tsv") {
                if (table_style == "granges") {
                  result_df <- as.data.frame(result_grl)
                  result_df <-
                    result_df[, colnames(result_df) != "group"]
                  colnames(result_df)[colnames(result_df) == "group_name"] = "mod_name"
                  result_df <-
                    cbind(result_df, as.data.frame(mcols(result_grl))[rep(seq_along(result_grl), elementNROWS(result_grl)),])

                } else {
                  scores <- -1 * log2(mcols(result_grl)$padj)

                  scores[is.na(scores)] <- 0

                  mcols(result_grl)$score <- scores


                  export(
                    object = result_grl,
                    con = file.path(save_dir, paste0(file_name, ".bed")),
                    format = "BED"
                  )

                  result_df <-
                    read.table(file.path(save_dir, paste0(file_name, ".bed")),
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

                  result_df$geneID <- sapply( result_grl, function(x) x$gene_id[1] )

                  mcols(result_grl) <-
                    mcols(result_grl)[,!colnames(mcols(result_grl)) %in% "score"]

                  result_df <-
                    cbind(result_df , as.data.frame(mcols(result_grl)))

                }


                write.table(
                  result_df,
                  file = file.path(save_dir, paste0(file_name, ".txt")),
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
                  con = file.path(save_dir, paste0(file_name, ".",  tolower(format))),
                  format = format
                )

              }
            }
          })
