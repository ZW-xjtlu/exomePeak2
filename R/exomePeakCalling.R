#' @title Exome peak calling on MeRIP-seq datasets while considering the biological variabilities.
#'
#' @description \code{exomePeakCalling} call RNA methylation peaks on exome regions with statistical tests that account for the biological variabilities between samples.
#'
#' @details The function conduct exome level peak calling based on the read alignment results and the transcript annotations.
#'
#' @param merip_bams a \code{MeripBamFileList} object.
#' @param txdb a \code{TxDb} object, it can also be a single character string such as "hg19" which is processed by \code{\link{makeTxDbFromUCSC}}.
#' @param bsgenome a \code{\link{BSgenome}} object for the genome sequence, alternatively it could be the name of the reference genome recognized by \code{\link{getBSgenom}}.
#' @param gene_anno a string, which specifies a gene annotation GFF/GTF file if available, default: NA.
#' @param fragment_length a positive integer of the expected fragment length in bp; default 100.
#' @param binding_length a positive integer of the antibody binding length in IP samples; default 25.
#' @param step_length a positive integer of the shift size of the sliding window; default is the binding length.
#' @param glm_type a character, which can be one of the "auto", "poisson", "NB", and "DESeq2". This argument specify the type of generalized linear model used in peak calling; Default to be "auto".
#'
#' Under the default setting, the DESeq2 method is implemented on experiments with more than 3 biological replicates for both IP and input samples.
#' The poisson GLM will be implemented otherwise.
#'
#' @param count_cutoff a non negative integer value of the minimum average reads count per window used in peak calling; default 5.
#' @param p_cutoff a value of the p value cut-off used in peak calling; default NULL.
#' @param p_adj_cutoff a value of the adjusted p value cutoff used in DESeq inference; default 0.05.
#' @param logFC_cutoff a non negative numeric value of the log2 fold change (log2 IP/input) cutoff used in the inferene of peaks.
#' @param peak_width positive integer of the minimum width for the merged peaks; default \code{fragment_length} .
#' @param drop_overlapped_genes a logical indicating whether the bins on overlapping genes are dropped or not; default TRUE.
#' @param parallel a logical indicating whether to use parallel computation, consider this if your computer has more than 16GB RAM.
#' @param mod_annotation a \code{GRanges} object for user provided single based RNA modification annotation. If provided, the peak calling step will be skipped.
#' @param mask_5p a logical of whether to exclude the transcript five prime regions in control; default TRUE.
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
                      txdb = NULL,
                      bsgenome = NULL,
                      gene_anno_gff = NULL,
                      fragment_length = 100,
                      binding_length = 25,
                      step_length = binding_length,
                      count_cutoff = 5,
                      p_cutoff = NULL,
                      p_adj_cutoff = 0.05,
                      logFC_cutoff = 0,
                      peak_width = fragment_length/2,
                      glm_type = c("auto","poisson","NB","DESeq2"),
                      drop_overlapped_genes = TRUE,
                      parallel = FALSE,
                      bp_param = NULL,
                      mod_annotation = NULL,
                      background = NULL,
                      mask_5p = TRUE
                      ){
  glm_type <- match.arg(glm_type)

  stopifnot(fragment_length > 0)

  stopifnot(step_length > 0)

  stopifnot(peak_width > 0)

  stopifnot(logFC_cutoff >= 0)

  stopifnot(count_cutoff >= 0)

  if(is.null(bsgenome)) {
    warning("bsgenome is not provided, peak calling will be performed without GC correction.\nIt is highly recommended to conduct GC conrrection for accurate peak calling and quantification.", call. = FALSE, immediate. = TRUE)
  }

  if(!is.null(gene_anno_gff)) {
    txdb <- makeTxDbFromGFF( gene_anno_gff )
  } else {
    if(is.null(txdb)) {
      stop("missing transcript annotation, please provide either txdb object or GFF/GTF file.")
    }

    if(!is(txdb,"TxDb")){
      txdb <- makeTxDbFromUCSC( txdb )
    }
  }

  paired <- any( asMates(merip_bams) )

  if(is.null(mod_annotation)) {

  message("generate bins on exons.")

  exome_bins_grl <- exome_bins_from_txdb(txdb = txdb,
                                        window_size = binding_length,
                                        step_size = step_length,
                                        drop_overlapped_genes = drop_overlapped_genes)

  message("count reads 5'TAGs on bins.")

  split_x <- function(x){return(split(x, names(x)))}

  if(!parallel) {
    register(SerialParam())
    register(MulticoreParam(workers = 1))
    register(SnowParam(workers = 1))
  } else {
    if(!is.null(bp_param)){
      register(bp_param,default = TRUE)
    }
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

  #calculate effective GC content if BSgenome is provided

  if(!is.null(bsgenome)){

   indx_gr <- names(unlist(rowRanges(SE_Peak_counts)))

   gc_freq <- as.vector( letterFrequency( Views(bsgenome, unlist(rowRanges(SE_Peak_counts))), letters="CG" ) )

   region_widths <- sum(width(rowRanges(SE_Peak_counts)))

   gc_contents <- tapply(gc_freq, indx_gr, sum) / region_widths

   rowData_SE <- DataFrame(gc_contents = gc_contents,
                        region_widths = region_widths)

   rm(indx_gr, gc_freq, region_widths, gc_contents)
  }

  rowRanges(SE_Peak_counts) <- exome_bins_grl[as.numeric(rownames(SE_Peak_counts))] #replace the row ranges with the not expanded version

  rm(exome_bins_grl) #release RAM

  #Replace the row ranges with only the binding regions, but the assays are the sum of the flanking regions.

  colData(SE_Peak_counts) = DataFrame(metadata(merip_bams))

  rowData(SE_Peak_counts) = rowData_SE

  rm(rowData_SE)

  #Peak calling with user defined statistical methods

  if(glm_type == "auto") {
    if(all( table(colData(SE_Peak_counts)$design_IP) > 3)) {
      glm_type <- "DESeq2"
    } else {
      glm_type <- "poisson"
    }
  }

  if(glm_type == "poisson") {
    message("peak calling with poisson GLM Wald test")
  }

  if(glm_type == "NB") {
    message("peak calling with NB GLM Wald test")
  }

  if(glm_type == "DESeq2") {
    message("peak calling with DESeq2 NB GLM Wald test")
  }

  #Directly return the merged summarized experiment
  gr_meth <- call_peaks_with_GLM( SE_bins = SE_Peak_counts,
                                  glm_type = glm_type,
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
                                  drop_overlapped_genes = drop_overlapped_genes,
                                  drop_5p = mask_5p,
                                  distance_5p = 200,
                                  control_width = peak_width)

  annotation_row_features <- count_row_features

  #Filter and rename the methylation peaks.

  indx_meth <- grepl("meth_", names( annotation_row_features ))

  annotation_row_features[indx_meth] <- gr_meth[ gsub("meth_", "", grep("meth_", names( count_row_features ) , value = TRUE) ) ]

  indx_keep <- !vector("logical", length = length(annotation_row_features))

  indx_keep[indx_meth][sum(width(annotation_row_features[indx_meth])) < peak_width] = FALSE

  annotation_row_features <- annotation_row_features[indx_keep]

  count_row_features <- count_row_features[indx_keep]

  rm(gr_meth, gr_meth_flanked, indx_meth, indx_keep)

  message("count reads on the merged peaks and the control regions")

  if(!parallel) {
    register(SerialParam(),default = FALSE)
    register(MulticoreParam(workers = 1))
    register(SnowParam(workers = 1))
  } else {
    if(!is.null(bp_param)){
      register(bp_param,default = TRUE)
    }
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
                                    drop_overlapped_genes = drop_overlapped_genes,
                                    drop_5p = mask_5p,
                                    distance_5p = 200,
                                    control_width = peak_width)

    message("count reads using single base annotation on exons")

    if(!parallel) {
      register(SerialParam())
      register(MulticoreParam(workers = 1))
      register(SnowParam(workers = 1))
    } else {
      if(!is.null(bp_param)){
        register(bp_param,default = TRUE)
      }
    }

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

  if(!is.null(bsgenome)){
    elementMetadata( SummarizedExomePeaks ) <- GC_content_over_grl(
      bsgenome = bsgenome,
      txdb = txdb,
      grl = rowRanges( SummarizedExomePeaks ),
      fragment_length = fragment_length,
      binding_length = binding_length,
      drop_overlapped_genes = drop_overlapped_genes,
      effective_GC = FALSE
    )
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
