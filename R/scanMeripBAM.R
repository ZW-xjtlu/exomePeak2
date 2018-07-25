#' @title Scan the BAM files of a MeRIP-seq experiment.
#'
#' @description \code{scanMeripBAM} is used to check and organize the BAM files in MeRIP-seq data before peak calling,
#' the flag parameters for the filtering can be defined at this step.
#'
#' @details The function takes the input of MeRIP-seq BAM files directories.
#' The provided BAM files are checked, and the design information of IP over input and treated over untreated are stored in the output.
#'
#' @param bam_ip a character vector of the BAM file directories for the IP samples.
#' @param bam_input a character vector of the BAM file directories for the input samples.
#' @param bam_treated_ip a character vector of the BAM file directories for the IP samples in treated group.
#' @param bam_treated_input a character vector of the BAM file directories for the input samples in treated group.
#'
#' If the dataset does not contain the extra design of treatment (such as the perturbation of the regulators), please fill only the \code{BAM_ip} and \code{BAM_input} arguments.
#'
#' @param paired_end a logical value indicating the library types, TRUE if the read is from paired end library, otherwise it will be treated as the single end reads.
#' @param random_primer a logical value indicating whether the library used is strand specific, default to be FALSE.
#' @param index_bam a logical value indicating whether to sort and index BAM automatically if the bam index is not found; default TRUE.
#'
#' The BAM index files will be named by adding ".bai" after the names of corresponding BAM files.
#'
#' @param bam_files an optional character string of all the BAM files to be analyzed, if it is provided, the first 4 arguments above will be ignored.
#' @param design_ip an optional logical vector indicating the design of IP and input, with TRUE represents for IP.
#' @param design_treatment an optional logical vector indicating the design of treatment and control, with TRUE represents for treated samples.
#' @param mapq A non-negative integer specifying the minimum mapping quality to include. BAM records with mapping qualities less than mapq are discarded; default 30L.
#' @param isSecondaryAlignment,isNotPassingQualityControls,isDuplicate,... arguments that determine the sam flag filtering, inherited from \code{\link{ScanBAMParam}}.
#'
#' @return This function will return a list with 2 named elements: BAMList and Parameter.
#'
#' @examples
#' MeRIP_Seq_Alignment <- scanMeripBAM(
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
#' MeRIP_Seq_Alignment2 <- scanMeripBAM(
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
#' @importFrom Rsamtools BamFileList ScanBamParam scanBamFlag path indexBam sortBam index<-
#' @export

scanMeripBAM <- function(bam_ip = NULL,
                     bam_input = NULL,
                     bam_treated_ip = NULL,
                     bam_treated_input = NULL,
                     paired_end = FALSE,
                     random_primer = TRUE,
                     index_bam = TRUE,
                     bam_files = NULL,
                     design_ip = NULL,
                     design_treatment = NULL,
                     mapq = 30L,
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


bai_temp = paste0(bam_files,".bai")

sorted_bai_temp = gsub(".bam$","_sorted.bam.bai",bam_files)

exist_indx <- all( file.exists( bai_temp ) ) | all( file.exists( sorted_bai_temp ) )

if(!exist_indx){

  if(!index_bam) {
    warning(paste0("cannot find the bam index files under: ",
                   paste0( index(bam.list)[!exist_indx], collapse = ", "),
                   ", The bam files are treated as not indexed."),
                   call. = F, immediate. = T )
  } else {
     message("The BAM files are not indexed, sort and indexing BAM files using RsamTools...")

     sorted_bam_names <- gsub( ".bam$", "_sorted", bam_files )

     for(i in seq_along(sorted_bam_names)){
      suppressWarnings( sortBam(bam_files[i], destination = sorted_bam_names[i]) )
     }

     indexBam(paste0( sorted_bam_names, ".bam" ))

     bam.list = BamFileList(
       file = paste0( sorted_bam_names , ".bam" ),
       asMates=paired_end
     )

    index(bam.list) = paste0(sorted_bam_names,".bam.bai")

  }

} else {
index(bam.list) = bai_temp
}

rm(bai_temp, sorted_bai_temp)


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

return(

new("MeripBamFileList",
    listData = bam.list@listData,
    elementType = bam.list@elementType,
    elementMetadata = bam.list@elementMetadata,
    metadata = bam.list@metadata,
    Parameter =  ScanBamParam(what = "mapq", flag = bam_flag, mapqFilter = mapq),
    RandomPrimer = random_primer)

)
}

