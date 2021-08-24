#' Method Parameter
#' @name Parameter-methods
#' @rdname Parameter-methods
#' @exportMethod Parameter
setGeneric("Parameter", function(x) {standardGeneric("Parameter")})

#' Method LibraryType
#' @name LibraryType-methods
#' @rdname LibraryType-methods
#' @exportMethod LibraryType
setGeneric("LibraryType", function(x) standardGeneric("LibraryType"))

#' Method exomePeakCalling
#' @name exomePeakCalling-methods
#' @rdname exomePeakCalling-methods
#' @exportMethod exomePeakCalling
setGeneric("exomePeakCalling", function(merip_bams = NULL,
                                        txdb = NULL,
                                        bsgenome = NULL,
                                        genome = NA,
                                        mod_annot = NULL,
                                        glm_type = c("DESeq2", "NB", "Poisson"),
                                        background_method = c("Gaussian_mixture",
                                                              "m6Aseq_prior",
                                                              "manual",
                                                              "all"),
                                        manual_background = NULL,
                                        correct_GC_bg = TRUE,
                                        qtnorm = FALSE,
                                        gff_dir = NULL,
                                        fragment_length = 100,
                                        binding_length = 25,
                                        step_length = binding_length,
                                        min_peak_width = fragment_length/2,
                                        max_peak_width = Inf,
                                        pc_count_cutoff = 5,
                                        bg_count_cutoff = 50,
                                        p_cutoff = 1e-05,
                                        p_adj_cutoff = NULL,
                                        log2FC_cutoff = 0,
                                        parallel = 3,
                                        bp_param = NULL) {standardGeneric("exomePeakCalling")})

#' Method estimateSeqDepth
#' @name estimateSeqDepth-methods
#' @rdname estimateSeqDepth-methods
#' @exportMethod estimateSeqDepth
setGeneric("estimateSeqDepth", function(sep,
                                        from = c("Background","Modification","All"),
                                        ...) {standardGeneric("estimateSeqDepth")})

#' Method plotSizeFactors
#' @name plotSizeFactors-methods
#' @rdname plotSizeFactors-methods
#' @exportMethod plotSizeFactors
setGeneric("plotSizeFactors", function(sep) {standardGeneric("plotSizeFactors")})

#' Method GCsizeFactors
#' @name GCsizeFactors-methods
#' @rdname GCsizeFactors-methods
#' @exportMethod GCsizeFactors
setGeneric("GCsizeFactors", function(x1) {standardGeneric("GCsizeFactors")})

#' Method GCsizeFactors<-
#' @name GCsizeFactors-methods
#' @rdname GCsizeFactors-methods
#' @exportMethod GCsizeFactors<-
setGeneric("GCsizeFactors<-", function(x2,value) {standardGeneric("GCsizeFactors<-")})

#' Method exomePeak2Results
#' @name exomePeak2Results-methods
#' @rdname exomePeak2Results-methods
#' @exportMethod exomePeak2Results
setGeneric("exomePeak2Results", function(x1) {standardGeneric("exomePeak2Results")})

#' Method exomePeak2Results<-
#' @name exomePeak2Results-methods
#' @rdname exomePeak2Results-methods
#' @exportMethod exomePeak2Results<-
setGeneric("exomePeak2Results<-", function(x2,value) {standardGeneric("exomePeak2Results<-")})

#' Method normalizeGC
#' @name normalizeGC-methods
#' @rdname normalizeGC-methods
#' @exportMethod normalizeGC
setGeneric("normalizeGC", function(sep,
                                   bsgenome = "hg19",
                                   txdb = "hg19",
                                   gff_dir = NULL,
                                   fragment_length = 100,
                                   binding_length = 25,
                                   feature = c("All","Background","Modification"),
                                   qtnorm = FALSE,
                                   effective_GC = FALSE) {standardGeneric("normalizeGC")})


#' Method glmM
#' @name glmM-methods
#' @rdname glmM-methods
#' @exportMethod glmM
setGeneric("glmM", function(sep,
                               glm_type = c("DESeq2","NB","Poisson"),
                               LFC_shrinkage = c("apeglm","Gaussian","ashr","none"),
                               ...) {standardGeneric("glmM")})

#' Method glmDM
#' @name glmDM-methods
#' @rdname glmDM-methods
#' @exportMethod glmDM
setGeneric("glmDM", function(sep,
                             glm_type = c("DESeq2","NB","Poisson"),
                             LFC_shrinkage = c("apeglm","ashr"),
                             ...) {standardGeneric("glmDM")})


# @rdname plotGuitar
# @export
# setGeneric("plotGuitar", function(sep,
#                                  txdb = NULL,
#                                  save_pdf_prefix = NULL,
#                                  include_control_regions = TRUE,
#                                  guitar_coordinate = NULL,
#                                  save_dir = ".") {standardGeneric("plotGuitar")})

#' Method plotExonLength
#' @name plotExonLength-methods
#' @rdname plotExonLength-methods
#' @exportMethod plotExonLength
setGeneric("plotExonLength", function(sep,
                                      txdb = NULL,
                                      save_pdf_prefix = NULL,
                                      include_control_regions = TRUE,
                                      save_dir = ".") {standardGeneric("plotExonLength")})

#' Method plotReadsGC
#' @name plotReadsGC-methods
#' @rdname plotReadsGC-methods
#' @exportMethod plotReadsGC
setGeneric("plotReadsGC", function(sep,
                                   bsgenome = NULL,
                                   txdb = NULL,
                                   save_pdf_prefix = NULL,
                                   fragment_length = 100,
                                   binding_length = 25,
                                   effective_GC = FALSE,
                                   pool_replicates = FALSE,
                                   save_dir = ".") {standardGeneric("plotReadsGC")})

#' Method plotLfcGC
#' @name plotLfcGC-methods
#' @rdname plotLfcGC-methods
#' @exportMethod plotLfcGC
setGeneric("plotLfcGC", function(sep,
                                  bsgenome = NULL,
                                  txdb = NULL,
                                  save_pdf_prefix = NULL,
                                  point_size = 0.05,
                                  xlim = c(0.2,0.9), 
                                  fragment_length = 100,
                                  binding_length = 25,
                                  effective_GC = FALSE,
                                  save_dir = ".") {standardGeneric("plotLfcGC")})

#' Method exportResults
#' @name exportResults-methods
#' @rdname exportResults-methods
#' @exportMethod exportResults
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

#' Method Results
#' @name Results-methods
#' @rdname Results-methods
#' @exportMethod Results
setGeneric("Results", function(sep,
                               cut_off_pvalue = NULL,
                               cut_off_padj = 0.05,
                               cut_off_log2FC = 0,
                               min_num_of_positive = 30,
                               expected_direction = c("both", "hyper", "hypo"),
                               inhibit_filter = FALSE,
                               table_style = c("bed","granges")) {standardGeneric("Results")})


