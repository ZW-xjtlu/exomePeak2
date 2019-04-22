#' @title Retrieve GC content level from genome over regions defined by a GRanges or GRangesList object.
#' @param grl A GRangesList object.
#' @param bsgenome A BSgenome object.
#' @param txdb A TxDb object.
#' @param fragment_length A positive integer of the expected fragment length in the RNA-Seq library; Default 100.
#' @param binding_length A positive integer of the expected antibody binding length of IP; Default 25.
#' @param effective_GC If TRUE, the GC content calculation will be weighted by the fragment mapping probabilities,
#' currently it is only supported for the single based modification annotation; Default FALSE.
#'
#' @description The GC content is calculated only on exon sequences, for regions that failed to mapped to exons in txdb, the function would return NA instead.
#' The ranges in grl will be flanked by the size fragment_length - binding_length before GC content calculation.
#'
#' @return a DataFrame contains 2 columns:
#'
#' 1. GC_content: the (effective) GC content of each GRanges / GRangesList element
#' 2. Indx_length: the calculated total feature length used to quantify reads.
#'
#' @import BSgenome
#' @import GenomicRanges
#' @importFrom S4Vectors DataFrame
#'
#'
GC_content_over_grl <- function(bsgenome,
                                txdb,
                                grl,
                                fragment_length = 100,
                                binding_length = 25,
                                effective_GC = FALSE) {

stopifnot(is(grl,"GRangesList"))

stopifnot(fragment_length > 0)

grl_mod <- grl[grepl("mod_",names(grl))]

all_sb <- all(sum(width(grl_mod)) == 1) & all(elementNROWS(grl_mod) == 1)

if(effective_GC & !all_sb) {
  stop("The effective GC content cannot be calculated on data that is not counted using single based annotation.")
}

if(all_sb) {
  message("Single based modification annotations detected, the flanking lengths are adjusted to floor( fragment_length - binding_length/2).")
}

if(!effective_GC){
#Calculate naive GC content.

flank_length <- ifelse( all_sb,
                        floor( fragment_length - binding_length / 2 ),
                        fragment_length - binding_length )

flanked_mod_gr <-  suppressWarnings( flank_on_exons( grl = grl_mod,
                                                     flank_length = flank_length,
                                                     txdb = txdb,
                                                     index_flank = FALSE ) )

flanked_mod_grl <- split( flanked_mod_gr, names( flanked_mod_gr ) )

rm(flanked_mod_gr)

flanked_grl <- c( flanked_mod_grl,
                  grl[ grepl("control_", names(grl)) ] )

rm(flanked_mod_grl)

flanked_gr <- unlist( flanked_grl )

names(flanked_gr) <- gsub("\\..*$", "", names(flanked_gr) )

GC_freq <- as.vector( letterFrequency(Views(bsgenome,flanked_gr), letters="CG") )

sum_freq <- tapply(GC_freq,names(flanked_gr), sum)

sum_freq <- sum_freq[names(flanked_grl)]

GC_return <- rep(NA, length(grl))

names(GC_return) <- names(grl)

GC_return[match(names(flanked_grl), names(grl))] <- sum_freq/sum(width(flanked_grl))

Width_return <- rep(NA, length(grl))

names(Width_return) <- names(grl)

Width_return[match(names(flanked_grl), names(grl))] <- sum(width(flanked_grl))

rm(flanked_gr,flanked_grl)

} else {

#Calculate weighted / effective GC content on genome coordinate.

flank <- fragment_length - binding_length

gcpos <- startIndex( vmatchPattern("S", DNAStringSet( Views(bsgenome, unlist(grl_mod) + (flank + floor(binding_length/2))) ), fixed = "subject") )

weight <- c(seq_len(flank), rep(flank + 1, binding_length), rev(seq_len(flank)))

weight <- weight/sum(weight)

GC_mod <- round(sapply(gcpos, function(x) sum(weight[x])), 3)

names(GC_mod) <- names(grl_mod)

#Calculate the background GC content using ordinary method without flanking.

GC_control <- as.vector( letterFrequency( Views(Hsapiens,unlist(grl[grepl("control_",names(grl))])) , letters = "GC", as.prob = T) )
names(GC_control) <- grep("control_",names(grl),value = T)

GC_return <- rep(NA,length(grl))

names(GC_return) <- names(grl)

GC_return[match(names(c(GC_mod,GC_control)), names(grl))] <- c(GC_mod,GC_control)

Width_return <- sum( width( grl ) )

Width_return[grepl("mod_", names(grl))] <- 1 + 2*floor( fragment_length - binding_length / 2 )

}

return(
  DataFrame(GC_content = GC_return,
            feature_length = Width_return)
       )

}
