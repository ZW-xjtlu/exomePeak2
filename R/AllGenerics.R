#' @rdname Parameter
#' @export
setGeneric("Parameter", function(x) {standardGeneric("Parameter")})

#' @rdname LibraryType
#' @export
setGeneric("LibraryType", function(x) standardGeneric("LibraryType"))

#' @rdname exomePeakCalling
#' @export
setGeneric("exomePeakCalling", function(merip_bams = NULL,
                                        txdb = NULL,
                                        bsgenome = NULL,
                                        mod_annot = NULL,
                                        glm_type = c("DESeq2", "NB", "Poisson"),
                                        background = c("Gaussian_mixture", "m6Aseq_prior", "manual", "all"),
                                        manual_background = NULL,
                                        correct_GC_bg = FALSE,
                                        qtnorm = TRUE,
                                        gff_dir = NULL,
                                        fragment_length = 100,
                                        binding_length = 25,
                                        step_length = binding_length,
                                        pc_count_cutoff = 5,
                                        bg_count_cutoff = 50,
                                        p_cutoff = NULL,
                                        p_adj_cutoff = 0.05,
                                        logFC_cutoff = 0,
                                        peak_width = fragment_length / 2,
                                        parallel = FALSE,
                                        bp_param = NULL) {standardGeneric("exomePeakCalling")})

#' @rdname estimateSeqDepth
#' @export
setGeneric("estimateSeqDepth", function(sep,
                                        from = c("Control","Modification","Both"),
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

#' @rdname normalizeGC
#' @export
setGeneric("normalizeGC", function(sep,
                                       bsgenome = "hg19",
                                       txdb = "hg19",
                                       gff_dir = NULL,
                                       fragment_length = 100,
                                       binding_length = 25,
                                       feature = c("background","all"),
                                       qtnorm = FALSE,
                                       effective_GC = FALSE) {standardGeneric("normalizeGC")})


#' @rdname glmM
#' @export
setGeneric("glmM", function(sep,
                               glm_type = c("auto","Poisson", "NB", "DESeq2"),
                               LFC_shrinkage = c("apeglm","Gaussian","ashr"),
                               ...) {standardGeneric("glmM")})

#' @rdname glmDM
#' @export
setGeneric("glmDM", function(sep,
                             glm_type = c("auto","Poisson", "NB", "DESeq2"),
                             LFC_shrinkage = c("apeglm","ashr"),
                             ...) {standardGeneric("glmDM")})


#' @rdname plotGuitar
#' @export
setGeneric("plotGuitar", function(sep,
                                  txdb = NULL,
                                  save_pdf_prefix = NULL,
                                  include_control_regions = TRUE,
                                  guitar_coordinate = NULL,
                                  save_dir = ".") {standardGeneric("plotGuitar")})

#' @rdname plotExonLength
#' @export
setGeneric("plotExonLength", function(sep,
                                      txdb = NULL,
                                      save_pdf_prefix = NULL,
                                      include_control_regions = TRUE,
                                      save_dir = ".") {standardGeneric("plotExonLength")})

#' @rdname plotReadsGC
#' @export
setGeneric("plotReadsGC", function(sep,
                                   bsgenome = NULL,
                                   txdb = NULL,
                                   save_pdf_prefix = NULL,
                                   fragment_length = 100,
                                   binding_length = 25,
                                   effective_GC = FALSE,
                                   pool_replicates = FALSE,
                                   save_dir = ".") {standardGeneric("plotReadsGC")})

#' @rdname plotLfcGC
#' @export
setGeneric("plotLfcGC", function(sep,
                                  bsgenome = NULL,
                                  txdb = NULL,
                                  save_pdf_prefix = NULL,
                                  fragment_length = 100,
                                  binding_length = 25,
                                  effective_GC = FALSE,
                                  save_dir = ".") {standardGeneric("plotLfcGC")})

#' @rdname exportResults
#' @export
setGeneric("exportResults", function(sep,
                                     format = c("CSV","BED","RDS"),
                                     table_style = c("bed","granges"),
                                     save_dir = "exomepeaks_result",
                                     cut_off_pvalue = NULL,
                                     cut_off_padj = 0.1,
                                     cut_off_log2FC = 0,
                                     min_num_of_positive = 30,
                                     expected_direction = c("both", "hyper", "hypo"),
                                     inhibit_filter = FALSE,
                                     reads_count = TRUE,
                                     GC_sizeFactors = TRUE) {standardGeneric("exportResults")})

#' @rdname Results
#' @export
setGeneric("Results", function(sep,
                               cut_off_pvalue = NULL,
                               cut_off_padj = 0.05,
                               cut_off_log2FC = 0,
                               min_num_of_positive = 30,
                               expected_direction = c("both", "hyper", "hypo"),
                               inhibit_filter = FALSE,
                               table_style = c("bed","granges")) {standardGeneric("Results")})


