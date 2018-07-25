#' @title Exome peak calling on MeRIP-seq datasets while considering the biological variabilities.
#'
#' @description \code{exomePeakCalling} call RNA methylation peaks on exome regions with statistical tests that account for the biological variabilities between samples.
#'
#' @details The function conduct exome level peak calling based on the read alignment results and the transcript annotations.
#'
#' @param merip_bams a \code{MeripBamFileList} object.
#' @param txdb a txdb object, it could be a single character string that can be recognized by \code{\link{makeTxDbFromUCSC}}; Default "hg19".
#' @param fragment_length a positive integer of the expected fragment length in bp; default 100.
#' @param binding_length a positive integer of the antibody binding length in IP samples; default 25.
#' @param step_length a positive integer of the shift size of the sliding window; default is the binding length.
#' @param count_cutoff a non negative integer value of the average reads count per window used in peak calling, default to be 5.
#' @param p_cutoff a value of the p value cut-off used in peak calling.
#' @param p_adj_cutoff a value of the adjusted p value cutoff used in DESeq inference; default 0.05.
#' @param logFC_cutoff a non negative numeric value of the log2 fold change (log2 IP/input) cutoff used in the inferene of peaks.
#' @param drop_overlapped_genes a logical indicating whether the bins on overlapping genes are dropped or not, default to be TRUE.
#' @param parallel a logical indicating whether to use parallel computation, consider this if your computer has more than 16GB RAM.
#' @param mod_annotation a \code{GRanges} object for user provided single based RNA modification annotation. If provided, the peak calling step will be skipped.
#' Reads count will be performed using the provided annotation flanked by length of floor(fragment_length - binding_length/2).
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
#' @importFrom GenomicAlignments summarizeOverlaps
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
                      txdb = "hg19",
                      fragment_length = 100,
                      binding_length = 25,
                      step_length = binding_length,
                      count_cutoff = 5,
                      p_cutoff = NULL,
                      p_adj_cutoff = 0.05,
                      logFC_cutoff = 0,
                      drop_overlapped_genes = TRUE,
                      parallel = FALSE,
                      mod_annotation = NULL,
                      background = NULL
                      ){

  stopifnot(fragment_length > 0)

  stopifnot(step_length > 0)

  stopifnot(logFC_cutoff >= 0)

  stopifnot(count_cutoff >= 0)

  if(is.null(txdb)) {
      stop("Transcript annotation undefined.")
  }

  if(!is(txdb,"TxDb")){
  txdb <- makeTxDbFromUCSC(txdb)
  }

  paired <- any( asMates(merip_bams) )

  if(is.null(mod_annotation)) {

  message("Extract bins on exons.")

  exome_bins_grl <- exome_bins_from_txdb(txdb = txdb,
                                        window_size = binding_length,
                                        step_size = step_length,
                                        drop_overlapped_genes = drop_overlapped_genes)

  message("Count reads on bins.")

  split_x <- function(x){return(split(x,names(x)))}

  if(!parallel) {
    register(SerialParam())
  }

  SE_Peak_counts <- summarizeOverlaps(
                    features = split_x(flank_on_exons(
                                       grl = exome_bins_grl,
                                       flank_length = fragment_length - binding_length,
                                       txdb = txdb,
                                       drop_overlapped_genes = drop_overlapped_genes,
                                       index_flank = FALSE
                                )),
                    reads = merip_bams,
                    param = Parameter(merip_bams),
                    mode = "Union",
                    inter.feature = FALSE,
                    preprocess.reads = reads_five_POS, #The reads are counted as the 5' POS
                    singleEnd = !paired,
                    ignore.strand = RandomPrimer(merip_bams),
                    fragments = paired
  )

  rowRanges(SE_Peak_counts) <- exome_bins_grl[as.numeric(rownames(SE_Peak_counts))] #replace the row ranges with the not expanded version

  rm(exome_bins_grl) #release RAM

  #Replace the row ranges with only the binding regions, but the assays are the sum of the flanking regions.

  colData(SE_Peak_counts) = DataFrame(metadata(merip_bams))

  #Peak calling with user defined statistical methods
  message("Peak calling with DESeq2 negative binomial Wald test.")

  #Directly return the merged summarized experiment
  gr_meth <- call_peaks_with_DESeq2( SE_bins = SE_Peak_counts,
                                     count_cutoff = count_cutoff,
                                     p_cutoff = p_cutoff,
                                     p_adj_cutoff = p_adj_cutoff,
                                     logFC_cutoff = logFC_cutoff,
                                     txdb = txdb,
                                     drop_overlapped_genes = drop_overlapped_genes )

  rm(SE_Peak_counts)

  #flank the merged methylation sites

  gr_meth_flanked <- flank_on_exons( grl = gr_meth,
                                     flank_length = fragment_length - binding_length,
                                     txdb = txdb,
                                     drop_overlapped_genes = drop_overlapped_genes,
                                     index_flank = FALSE )

  #calculate control regions that are disjoint of the methylation peaks

  count_row_features <- annot_bg( annot = gr_meth_flanked,
                                  txdb = txdb,
                                  cut_off_width = 1e5,
                                  cut_off_num = 2000,
                                  drop_overlapped_genes = drop_overlapped_genes )

  annotation_row_features <- count_row_features

  annotation_row_features[grepl("meth_", names( annotation_row_features ))] <- gr_meth[ gsub("meth_", "", grep("meth_", names( count_row_features ) , value = TRUE) ) ]

  rm(gr_meth, gr_meth_flanked)

  message("Count reads on the merged peaks and the control regions.")

  if(!parallel) {
    register(SerialParam())
  }

  SummarizedExomePeaks <- summarizeOverlaps(
                          features = count_row_features,
                          reads = merip_bams,
                          param = Parameter(merip_bams),
                          mode = "Union",
                          inter.feature = FALSE,
                          preprocess.reads = reads_five_POS, #The reads are counted as the 5' POS
                          singleEnd = !paired,
                          ignore.strand = RandomPrimer(merip_bams),
                          fragments = paired
                       )

  rowRanges(SummarizedExomePeaks) <- annotation_row_features

  rm(count_row_features, annotation_row_features)

  colData(SummarizedExomePeaks) <- DataFrame(metadata(merip_bams))

  } else {

    stopifnot( is( mod_annotation, "GRanges") | is( mod_annotation, "GRangesList") )

    stopifnot( all(width(mod_annotation)) == 1 )

    if(

      is( mod_annotation, "GRanges")

    ) {

      mod_annotation <- split( mod_annotation , seq_along(mod_annotation) )

    } else {

      names(mod_annotation) <- seq_along(mod_annotation) #make sure the annotation is indexed by integer sequence

    }

      mod_annotation_flanked <- flank_on_exons(grl = mod_annotation,
                                               flank_length = floor( fragment_length - binding_length/2 ),
                                               txdb = txdb,
                                               drop_overlapped_genes = drop_overlapped_genes,
                                               index_flank = FALSE)

      merged_peaks_grl <- annot_bg( annot = mod_annotation_flanked,
                                    txdb = txdb,
                                    cut_off_width = 1e5,
                                    cut_off_num = 2000,
                                    drop_overlapped_genes = drop_overlapped_genes )

    if(!parallel){
      register(SerialParam())
    }

    message("Count reads using single base annotation on exons.")

    SE_temp <- summarizeOverlaps(
      features = merged_peaks_grl,
      reads = merip_bams,
      param = Parameter(merip_bams),
      mode = "Union",
      inter.feature = FALSE,
      preprocess.reads = reads_five_POS,
      singleEnd = !paired,
      ignore.strand = RandomPrimer(merip_bams),
      fragments = paired
    )

    #Replace the rowRanges with the single based GRangesList.
    index_meth <- grepl("meth_",rownames(SE_temp))

    index_sb_annot <- as.numeric( gsub("meth_", "", rownames(SE_temp)[index_meth]) )

    meth_count <- assay(SE_temp)[index_meth,][match(as.numeric(names(mod_annotation)),index_sb_annot),]

    rownames(meth_count) <- paste0("meth_", seq_len(nrow(meth_count)))

    #Fill the rows that are not on exons with zero counts.
    meth_count[is.na(meth_count)] <- 0

    control_count <- assay(SE_temp)[!index_meth,]

    #replace the metadata collumn with gene_ids
    mod_annotation_gr <- unlist(mod_annotation)

    mcols(mod_annotation_gr) <- DataFrame(gene_id = NA)

    mcols(mod_annotation_gr)$gene_id[index_sb_annot] <- mcols(unlist(rowRanges(SE_temp)))$gene_id[index_sb_annot]

    mod_annotation <- split(mod_annotation_gr,names(mod_annotation_gr))

    mod_annotation <- mod_annotation[order(as.numeric(names(mod_annotation)))]

    names(mod_annotation) <- paste0("meth_",names(mod_annotation))

    SummarizedExomePeaks <- SummarizedExperiment(

      assay = rbind(meth_count,control_count),

      rowRanges = c(mod_annotation,
                    rowRanges(SE_temp)[!index_meth]),

      colData = DataFrame(metadata(merip_bams))

    )

  }

  if(!is.null(background)) {

    rowRanges(SummarizedExomePeaks) <- replace_bg( grl = rowRanges(SummarizedExomePeaks),
                                                   bg = background,
                                                   txdb = txdb,
                                                   drop_overlapped_genes = drop_overlapped_genes )
  }

  return(

    new("SummarizedExomePeak",
        rowRanges = SummarizedExomePeaks@rowRanges,
        colData = SummarizedExomePeaks@colData,
        assays = SummarizedExomePeaks@assays,
        NAMES = SummarizedExomePeaks@NAMES,
        elementMetadata = SummarizedExomePeaks@elementMetadata,
        metadata = SummarizedExomePeaks@metadata,
        DESeq2Results = data.frame() )
  )

}

)
