#' @title Organize the BAM Files Information of a MeRIP-seq Data Set.
#'
#' @description \code{scanMeripBAM} check and organize the BAM files in MeRIP-seq data before peak calling using \code{\link{exomePeakCalling}}.
#' The library types of the RNA-seq and the filters such as SAM FLAG score are specified in this function.
#'
#' @details \code{scanMeripBAM} takes the input of the BAM file directories for the MeRIP-seq datasets.
#' It first checks the completeness of the BAM files and the BAM indexes. Then, the design information of IP/input and treated/control are returned as a \code{MeripBamFileList} object.
#' If the BAM file indexes are missing, the BAM files will be automatically indexed with the package \code{Rsamtools}.
#'
#' @param bam_ip a \code{character} vector for the BAM file directories of the (control) IP samples.
#' @param bam_input a \code{character} vector for the BAM file directories of the (control) input samples.
#' @param bam_treated_ip a \code{character} vector for the BAM file directories of the treated IP samples.
#' @param bam_treated_input a \code{character} vector for the BAM file directories of the treated input samples.
#'
#' If the bam files do not contain treatment group, user should only fill the arguments of \code{BAM_ip} and \code{BAM_input}.
#'
#' @param paired_end a \code{logical} of whether the data comes from the Paired-End Library, \code{TRUE} if the data is Paired-End sequencing; default \code{= FALSE}.
#' @param library_type a \code{character} specifying the protocal type of the RNA-seq library, can be one in \code{c("unstranded", "1st_strand", "2nd_strand")}; default \code{= "unstranded"}.
#'
#' \describe{
#' \item{\strong{unstranded}}{The randomly primed RNA-seq library type, i.e. both the strands generated during the first and the second strand sythesis are sequenced; example: Standard Illumina.}
#' \item{\strong{1st_strand}}{The first strand-specific RNA-seq library, only the strand generated during the first strand sythesis is sequenced; examples: dUTP, NSR, NNSR.}
#' \item{\strong{2nd_strand}}{The second strand-specific RNA-seq library, only the strand generated during the second strand sythesis is sequenced; examples: Ligation, Standard SOLiD.}
#' }
#'
#' @param index_bam a \code{logical} value indicating whether to sort and index BAM files automatically if the bam indexes are not found; default \code{= TRUE}.
#'
#' The BAM index files will be named by adding ".bai" after the names of the corresponding BAM files.
#'
#' @param bam_files optional, a \code{character} vector for all the BAM file directories, if it is provided, the first 4 arguments above will be ignored.
#'
#' @param design_ip optional, a \code{logical} vector indicating the information of IP/input, \code{TRUE} represents IP samples.
#'
#' @param design_treatment optional, a \code{logical} vector indicating the design of treatment/control, \code{TRUE} represents treated samples.
#'
#' @param mapq a non-negative integer specifying the minimum reads mapping quality. BAM records with mapping qualities less than \code{mapq} are discarded; default \code{= 30L}.
#'
#' @param isSecondaryAlignment,isNotPassingQualityControls,isDuplicate,... arguments specifying the filters on SAM FLAG scores, inherited from \code{\link{ScanBAMParam}}.
#'
#' @return a \code{MeripBamFileList} object.
#'
#' @examples
#'
#' ###For MeRIP-seq experiment without treatment group
#'
#' MeRIP_Seq_Alignment <- scanMeripBAM(
#'   bam_ip = c("IP_rep1.bam",
#'              "IP_rep2.bam",
#'              "IP_rep3.bam"),
#'   bam_input = c("input_rep1.bam",
#'                 "input_rep2.bam",
#'                 "input_rep3.bam"),
#'   paired_end = TRUE
#' )
#'
#' ###For MeRIP-seq experiment with treatment group
#'
#' MeRIP_Seq_Alignment <- scanMeripBAM(
#' bam_ip = c("IP_rep1.bam",
#'            "IP_rep2.bam",
#'            "IP_rep3.bam"),
#' bam_input = c("input_rep1.bam",
#'               "input_rep2.bam",
#'               "input_rep3.bam"),
#' bam_treated_ip = c("IP_treated_rep1.bam",
#'                    "IP_treated_rep2.bam"),
#' bam_treated_input = c("input_treated_rep1.bam",
#'                       "input_treated_rep2.bam"),
#' paired_end = TRUE
#' )
#'
#' @seealso \code{\link{exomePeakCalling}}
#' @import Rsamtools
#' @export

scanMeripBAM <- function(bam_ip = NULL,
                         bam_input = NULL,
                         bam_treated_ip = NULL,
                         bam_treated_input = NULL,
                         paired_end = FALSE,
                         library_type = c("unstranded","1st_strand","2nd_strand"),
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

  library_type <- match.arg(library_type)

  #Create bamfile list
  if (is.null(bam_files)) {
    bam_files <- c(bam_ip,
                   bam_input,
                   bam_treated_ip,
                   bam_treated_input)
  }

  exist_bam <- file.exists(bam_files)

  if( any(!exist_bam) ){
    stop(paste0("Files do not exist for:\n",
          paste(bam_files[!exist_bam],collapse = ", "),
         "\nplease check the input directories of the BAM files."))
  }

  rm(exist_bam)

  bam.list = BamFileList(file = bam_files,
                         asMates = paired_end)

  bai_temp = paste0(bam_files, ".bai")

  sorted_bai_temp = gsub(".bam$", "_sorted.bam.bai", bam_files)

  exist_indx <-
    all(file.exists(bai_temp)) |
    (all(file.exists(sorted_bai_temp)) &
       all(file.exists(
         gsub(".bam$", "_sorted.bam", bam_files)
       )))

  if (!exist_indx) {
    if (!index_bam) {
      warning(
        paste0(
          "Cannot find the bam index files, The bam files are treated as not indexed."
        ),
        call. = F,
        immediate. = T
      )
    } else {
      message("Cannot find the bam index files, sorting and indexing BAM files with Rsamtools...")

      sorted_bam_names <- gsub(".bam$", "_sorted", bam_files)

      for (i in seq_along(sorted_bam_names)) {
        suppressWarnings(sortBam(bam_files[i], destination = sorted_bam_names[i]))
      }

      indexBam(paste0(sorted_bam_names, ".bam"))

      bam.list = BamFileList(file = paste0(sorted_bam_names, ".bam"),
                             asMates = paired_end)

      index(bam.list) = normalizePath(paste0(sorted_bam_names, ".bam.bai"))

    }

  } else {
    if (all(file.exists(sorted_bai_temp))) {
      bam.list = BamFileList(file = gsub(".bam$", "_sorted.bam", bam_files),
                             asMates = paired_end)

      index(bam.list) = normalizePath(sorted_bai_temp)

    } else {
      index(bam.list) = normalizePath(bai_temp)

    }
  }

  rm(bai_temp, sorted_bai_temp)


  #Check the existence of the bam files
  exist_indx <- file.exists(path(bam.list))

  if (any(!exist_indx)) {
    stop(paste0("cannot find bam files under: ",
                paste0(path(bam.list)[!exist_indx] , collapse = ", ")))
  }

  #Create metadata of the bam files for experimental design information
  if (is.null(design_ip)) {
    design_ip = vector(length = length(bam_files))
    names(design_ip) = bam_files
    design_ip[c(bam_ip, bam_treated_ip)] = T
    names(design_ip) = NULL
  }

  if (is.null(design_treatment)) {
    design_treatment = vector(length = length(bam_files))
    names(design_treatment) = bam_files
    design_treatment[c(bam_treated_ip, bam_treated_input)] = T
    names(design_treatment) = NULL
  }

  stopifnot(length(design_ip) == length(design_treatment))

  mdf = data.frame(design_IP = design_ip,
                   design_Treatment = design_treatment)

  metadata(bam.list) <- mdf

  #Check if there are any duplicated bam names.
  dup_indx <- duplicated(names(bam.list))

  if (any(dup_indx)) {
    warning("Containing duplicated bam file names, the bam files are re-named.")
    ip_char <- rep("input", length(bam_files))
    ip_char[mdf$design_IP] <- "IP"
    trt_char <- rep("ctrl", length(bam_files))
    trt_char[mdf$design_Treatment] <- "Trt"
    new_name <- paste0(ip_char, "_", trt_char)
    names(bam.list) <-
      paste0(new_name, "_rep", sequence(as.integer(table(new_name)[unique(new_name)])))
  }

  #Finally create the attribute for bam flag

  bam_flag <- scanBamFlag(
    isPaired = isPaired,
    isProperPair = isProperPair,
    hasUnmappedMate = hasUnmappedMate,
    isSecondaryAlignment = isSecondaryAlignment,
    isNotPassingQualityControls = isNotPassingQualityControls,
    isDuplicate = isDuplicate,
    ...
  )

  return(
    new(
      "MeripBamFileList",
      listData = bam.list@listData,
      elementType = bam.list@elementType,
      elementMetadata = bam.list@elementMetadata,
      metadata = bam.list@metadata,
      Parameter =  ScanBamParam(
        what = "mapq",
        flag = bam_flag,
        mapqFilter = mapq
      ),
      LibraryType = library_type
    )
  )
}

