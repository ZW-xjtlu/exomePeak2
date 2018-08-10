#' @rdname Parameter
#' @export
setGeneric("Parameter", function(x) {standardGeneric("Parameter")})

#' @rdname RandomPrimer
#' @export
setGeneric("RandomPrimer", function(x) standardGeneric("RandomPrimer"))

#' @rdname exomePeakCalling
#' @export
setGeneric("exomePeakCalling", function(merip_bams = NULL,
                                        txdb = NULL,
                                        bsgenome = NULL,
                                        gene_anno_gff = NULL,
                                        fragment_length = 100,
                                        binding_length = 25,
                                        step_length = binding_length,
                                        glm_type = c("auto","poisson","NB","DESeq2"),
                                        count_cutoff = 5,
                                        p_cutoff = NULL,
                                        p_adj_cutoff = 0.05,
                                        logFC_cutoff = 0,
                                        peak_width = fragment_length/2,
                                        drop_overlapped_genes = TRUE,
                                        parallel = FALSE,
                                        mod_annotation = NULL,
                                        background = NULL,
                                        mask_5p = TRUE) {standardGeneric("exomePeakCalling")})

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
setGeneric("normalizeGC", function(sep,
                                       bsgenome = "hg19",
                                       txdb = "hg19",
                                       gene_anno_gff = NULL,
                                       fragment_length = 100,
                                       binding_length = 25,
                                       feature = c("background","all"),
                                       qtnorm = FALSE,
                                       effective_GC = FALSE,
                                       drop_overlapped_genes = TRUE) {standardGeneric("normalizeGC")})


#' @rdname glmMeth
#' @export
setGeneric("glmMeth", function(sep,
                               glm_type = c("auto","poisson", "NB", "DESeq2"),
                               shrinkage_method = c("apeglm","normal","ashr"),
                               ...) {standardGeneric("glmMeth")})

#' @rdname glmDM
#' @export
setGeneric("glmDM", function(sep,
                             glm_type = c("auto","poisson", "NB", "DESeq2"),
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
                                   bsgenome = NULL,
                                   txdb = NULL,
                                   save_pdf_prefix = NULL,
                                   fragment_length = 100,
                                   binding_length = 25,
                                   effective_GC = FALSE,
                                   drop_overlapped_genes = TRUE,
                                   pool_replicates = FALSE) {standardGeneric("plotReadsGC")})

#' @rdname plotBetaGC
#' @export
setGeneric("plotEffectGC", function(sep,
                                  bsgenome = NULL,
                                  txdb = NULL,
                                  save_pdf_prefix = NULL,
                                  fragment_length = 100,
                                  binding_length = 25,
                                  effective_GC = FALSE,
                                  drop_overlapped_genes = TRUE) {standardGeneric("plotEffectGC")})

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


