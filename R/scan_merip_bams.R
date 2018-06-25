#' @title Scan the bam files of a MeRIP-seq experiment.
#'
#' @description \code{scan_merip_bams} is used to check and organize the bam files in MeRIP-seq data before peak calling,
#' the flag parameters for the filtering can be defined at this step.
#'
#' @details The function takes the input of MeRIP-seq bam files directories.
#' The provided bam files are checked, and the design information of IP over input and treated over untreated are stored in the output.
#'
#' @param bam_ip a character vector of the bam file directories for the IP samples.
#' @param bam_input a character vector of the bam file directories for the input samples.
#' @param bam_treated_ip a character vector of the bam file directories for the IP samples in treated group.
#' @param bam_treated_input a character vector of the bam file directories for the input samples in treated group.
#'
#' If the dataset does not contain the extra design of treatment (such as the perturbation of the regulators), please fill only the \code{bam_ip} and \code{bam_input} arguments.
#'
#' @param paired_end a logical value indicating the library types, TRUE if the read is from paired end library, otherwise it will be treated as the single end reads.
#' @param strand_specific a logical value indicating whether the library used is strand specific, default to be FALSE.
#' @param bam_index a logical value indicating the existence of bam indexes, TRUE if the bam files are indexed.
#'
#' The bam index files should be named by adding ".bai" after the names of corresponding bam files.
#' If the index files cannot be found in the provided directory, the bam files would be loaded as if they are not indexed.
#'
#' @param bam_files an optional character string of all the bam files to be analyzed, if it is provided, the first 4 arguments above will be ignored.
#' @param design_ip an optional logical vector indicating the design of IP and input, with TRUE represents for IP.
#' @param design_treatment an optional logical vector indicating the design of treatment and control, with TRUE represents for treated samples.
#' @param isSecondaryAlignment,isNotPassingQualityControls,isDuplicate,... arguments that determine the sam flag filtering, inherited from \code{\link{ScanBamParam}}.
#'
#' @return This function will return a list with 2 named elements: BamList and Parameter.
#'
#' @examples
#' MeRIP_Seq_Alignment <- scan_merip_bams(
#'                             bam_ip = c("./bam/SRR1182619.bam",
#'                                      "./bam/SRR1182621.bam",
#'                                      "./bam/SRR1182623.bam"),
#'                             bam_input = c("./bam/SRR1182620.bam",
#'                                           "./bam/SRR1182622.bam",
#'                                           "./bam/SRR1182624.bam"),
#'                             bam_treated_ip = c("./bam/SRR1182603.bam",
#'                                                "./bam/SRR1182605.bam"),
#'                             bam_treated_input = c("./bam/SRR1182604.bam",
#'                                                   "./bam/SRR1182606.bam"),
#'                             paired_end = TRUE
#'                             )
#'
#' ###It will provide identical result with the following arguments:
#'
#' MeRIP_Seq_Alignment2 <- scan_merip_bams(
#'                              bam_files = c("./bam/SRR1182619.bam",
#'                                            "./bam/SRR1182621.bam",
#'                                            "./bam/SRR1182623.bam",
#'                                            "./bam/SRR1182620.bam",
#'                                            "./bam/SRR1182622.bam",
#'                                            "./bam/SRR1182624.bam",
#'                                            "./bam/SRR1182603.bam",
#'                                            "./bam/SRR1182605.bam",
#'                                            "./bam/SRR1182604.bam",
#'                                            "./bam/SRR1182606.bam"),
#'                              design_ip = metadata(MeRIP_Seq_Alignment[[1]])$design_IP,
#'                              design_treatment = metadata(MeRIP_Seq_Alignment[[1]])$design_Treatment,
#'                              paired_end = TRUE
#'                              )
#'
#'identical(MeRIP_Seq_Alignment, MeRIP_Seq_Alignment2)
#'
#' @seealso \code{\link{merip_peak_calling}}
#' @importFrom Rsamtools BamFileList ScanBamParam scanBamFlag path
#' @export

scan_merip_bams <- function(bam_ip = NULL,
                     bam_input = NULL,
                     bam_treated_ip = NULL,
                     bam_treated_input = NULL,
                     paired_end = FALSE,
                     strand_specific = FALSE,
                     bam_index = FALSE,
                     bam_files = NULL,
                     design_ip = NULL,
                     design_treatment = NULL,
                     isSecondaryAlignment = FALSE,
                     isNotPassingQualityControls = FALSE,
                     isDuplicate = FALSE,
                     isPaired = NA,
                     isProperPair = NA,
                     hasUnmappedMate = NA,
                     ...) {

#Create bamfile list
if(is.null(bam_files)) {
bam_files <- c(bam_ip,
               bam_input,
               bam_treated_ip,
               bam_treated_input)
}

bam.list = BamFileList(
                       file = bam_files,
                       asMates=paired_end
          )

if(bam_index) {
bai_temp = paste0(bam_files,".bai")

exist_indx <- file.exists( bai_temp )

if(any(!exist_indx)){
  warning(paste0("cannot find the bam index files under: ",
              paste0( index(bam.list)[!exist_indx] ,collapse = ", "),
              ", The bam files are treated as not indexed."),
          call. = F,immediate. = T)
} else {
index(bam.list) = bai_temp
}

}

#Check the existence of the bam files
exist_indx <- file.exists( path(bam.list) )

if(any(!exist_indx)){
  stop(paste0("cannot find bam files under: ",
              paste0( path(bam.list)[!exist_indx] ,collapse = ", ")))
}

#Create metadata of the bam files for experimental design information
if(is.null(design_ip)) {
  design_ip = vector(length = length(bam_files))
  names(design_ip) = bam_files
  design_ip[c(bam_ip,bam_treated_ip)] = T
  names(design_ip) = NULL
}

if(is.null(design_treatment)) {
  design_treatment = vector(length = length(bam_files))
  names(design_treatment) = bam_files
  design_treatment[c(bam_treated_ip,bam_treated_input)] = T
  names(design_treatment) = NULL
}

stopifnot(length(design_ip) == length(design_treatment))

mdf = data.frame(design_IP = design_ip,
                 design_Treatment = design_treatment)

metadata(bam.list) <- mdf

#Check if there are any duplicated bam names.
dup_indx <- duplicated(names(bam.list))

if(any(dup_indx)) {
  warning("Containing duplicated bam file names, the bam files are re-named.")
  ip_char <- rep("input",length(bam_files))
  ip_char[mdf$design_IP] <- "IP"
  trt_char <- rep("ctrl",length(bam_files))
  trt_char[mdf$design_Treatment] <- "Trt"
  new_name <- paste0(ip_char,"_",trt_char)
  names(bam.list) <- paste0(new_name,"_rep", sequence( as.integer( table(new_name)[unique(new_name)] ) ) )
}

#Finally create the attribute for bam flag

bam_flag <- scanBamFlag(isPaired = isPaired,
isProperPair = isProperPair,
hasUnmappedMate = hasUnmappedMate,
isSecondaryAlignment = isSecondaryAlignment,
isNotPassingQualityControls = isNotPassingQualityControls,
isDuplicate = isDuplicate,
...)

merip_bam_list <- list(BamList = bam.list,
     Parameter = ScanBamParam(flag = bam_flag),
     StrandSpecific = strand_specific)

return(merip_bam_list)
}

