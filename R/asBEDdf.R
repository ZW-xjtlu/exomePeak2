#'@title coerce granges list into bed like dataframe
#'@importFrom rtracklayer asBED
#'@import GenomicRanges
asBEDdf <- function(grl){
mcols_df <- as.data.frame( mcols(grl) )
scores <- grl$score
mcols(grl) <- NULL
gr_bed <- asBED(grl)
df_bed <- data.frame(chr = seqnames( gr_bed ),
                     chromStart = start(gr_bed),
                     chromEnd = end(gr_bed),
                     name = gr_bed$name,
                     score = scores,
                     strand = strand( gr_bed ),
                     thickStart = start(gr_bed),
                     thickEnd = end(gr_bed),
                     itemRgb = 0,
                     blockCount = elementNROWS( gr_bed$blocks ),
                     blockSizes = width(gr_bed$blocks),
                     blockStarts = start(gr_bed$blocks))
}
