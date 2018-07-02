#' @rdname Parameter
#' @export
setGeneric("Parameter", function(x) {standardGeneric("Parameter")})

#' @rdname RandomPrimer
#' @export
setGeneric("RandomPrimer", function(x) standardGeneric("RandomPrimer"))

#' @rdname exomePeakCalling
#' @export
setGeneric("exomePeakCalling",function(merip_bams = NULL,
                                       genome = "hg19",
                                       gff_dir = NULL,
                                       txdb = NULL,
                                       window_size = 100,
                                       step_size = 10,
                                       count_cutoff = 10,
                                       p_cutoff = NULL,
                                       p_adj_cutoff = 0.05,
                                       logFC_cutoff = 0,
                                       drop_overlapped_genes = TRUE,
                                       parallel = FALSE,
                                       mod_annotation = NULL,
                                       background = NULL) {standardGeneric("exomePeakCalling")})

#' @rdname estimateSeqDepth
#' @export
setGeneric("estimateSeqDepth", function(sep,
                                           from = c("Control","Methylation","Both"),
                                           ...) {standardGeneric("estimateSeqDepth")})

#' @rdname plotSizeFactors
#' @export
setGeneric("plotSizeFactors", function(sep) {standardGeneric("plotSizeFactors")})

#' @rdname GCsizeFactors
#' @export
setGeneric("GCsizeFactors", function(x) {standardGeneric("GCsizeFactors")})

#' @rdname GCsizeFactors<-
#' @export
setGeneric("GCsizeFactors<-", function(x,...,value) {standardGeneric("GCsizeFactors<-")})

#' @rdname DESeq2Results
#' @export
setGeneric("DESeq2Results", function(x) {standardGeneric("DESeq2Results")})

#' @rdname DESeq2Results<-
#' @export
setGeneric("DESeq2Results<-", function(x,...,value) {standardGeneric("DESeq2Results<-")})

#' @rdname GCnormalization
#' @export
setGeneric("GCnormalization", function(sep,
                              bsgenome = "hg19",
                              feature = c("background","all"),
                              qtnorm = FALSE,
                              fragment_length = 100,
                              glm_offset = TRUE) {standardGeneric("GCnormalization")})


#' @rdname glmMeth
#' @export
setGeneric("glmMeth", function(sep,
                               shrinkage_method = c("apeglm","normal","ashr"),
                               ...) {standardGeneric("glmMeth")})

#' @rdname glmDM
#' @export
setGeneric("glmDM", function(sep,
                             shrinkage_method = c("apeglm","ashr"),
                             ...) {standardGeneric("glmDM")})


#' @rdname plotGuitar
#' @export
setGeneric("plotGuitar", function(sep,
                                  txdb = NULL,
                                  save_pdf_prefix = NULL,
                                  include_control_regions = TRUE,
                                  guitar_coordinate = NULL) {standardGeneric("plotGuitar")})

#' @rdname plotExonLength
#' @export
setGeneric("plotExonLength", function(sep,
                                      txdb = NULL,
                                      save_pdf_prefix = NULL,
                                      include_control_regions = TRUE) {standardGeneric("plotExonLength")})

#' @rdname plotReadsGC
#' @export
setGeneric("plotReadsGC", function(sep,
                                   bsgenome,
                                   fragment_length = 100,
                                   save_pdf_prefix = NULL,
                                   combine_replicates = FALSE) {standardGeneric("plotReadsGC")})

#' @rdname plotBetaGC
#' @export
setGeneric("plotBetaGC", function(sep,
                                  bsgenome,
                                  save_pdf_prefix = NULL,
                                  fragment_length = 100) {standardGeneric("plotBetaGC")})

#' @rdname exportResults
#' @export
setGeneric("exportResults", function(sep,
                                     format = c("txt","BED","RDS"),
                                     file_name = "exomepeaks_result",
                                     cut_off_pvalue = NULL,
                                     cut_off_padj = 0.05,
                                     cut_off_log2FC = 0,
                                     min_num_of_positive = 2000,
                                     expected_direction = "both",
                                     inhibit_filter = FALSE,
                                     table_style = c("bed","granges")) {standardGeneric("exportResults")})

#' @rdname Results
#' @export
setGeneric("Results", function(sep,
                               cut_off_pvalue = NULL,
                               cut_off_padj = 0.05,
                               cut_off_log2FC = 0,
                               min_num_of_positive = 2000,
                               expected_direction = "both",
                               inhibit_filter = FALSE,
                               table_style = c("bed","granges")) {standardGeneric("Results")})


