#' @title Retrieve GC content level from genome over regions defined by a GRanges or GRangesList object.
#' @param grl a GRanges or GRangesList object.
#' @param bsgenome a BSgenome object.
#' @param fragment_length an integer of the fragment length.
#' @description if the GRanges contain regions less than fragment_length, it will be re-sized to fragment_length by fixing at center.
#' This is currently not supported for user provided GRangesList annotations.
#' @return a list contains 2 elements:
#'
#' 1. GC_content: the GC content of each GRanges / GRangesList element
#' 2. Indx_length: the calculated total feature length used to quantify reads.
#'
#' @import BSgenome
#' @import GenomicRanges
#' @export
GC_content_over_grl <- function(grl,
                                bsgenome,
                                fragment_length = 100){

stopifnot(is(grl,"GRangesList")|is(grl,"GRanges"))

if(is(grl,"GRangesList")){

gr <- unlist(grl)
GC_freq <- letterFrequency(Views(bsgenome,gr),letters="CG")
width_each_group <- sum(width(grl))
sum_freq <- tapply(GC_freq,names(gr),sum)
GC_return <- sum_freq/width_each_group

} else {

index_resize <- width(grl) < fragment_length
grl[index_resize] <- resize(grl[index_resize],fragment_length,fix = "center")
GC_return <- letterFrequency(Views(bsgenome,grl),letters="CG",as.prob = TRUE)
width_each_group <- width(grl)

}

return(
  list(GC_content = GC_return,
       Indx_length = width_each_group)
       )
}
