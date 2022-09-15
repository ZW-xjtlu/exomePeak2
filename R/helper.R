## A function to extract exons grouped by unambiguous genes
exonsByiGenes <- function(txdb){
  exbg <- exonsBy(txdb, by = "gene")
  exbg <- exbg[elementNROWS(range(exbg)) == 1]
  fol <- findOverlaps(exbg)
  fol <- fol[queryHits(fol) != subjectHits(fol)]
  ol_indx_M <- as.matrix(fol)
  if (nrow(ol_indx_M) == 0) {
    return(reduce(exbg))
  }
  else {
    rm(fol)
    new_gene_names_temp <- names(exbg)
    new_gene_names_list <- split(new_gene_names_temp, seq_along(new_gene_names_temp))
    for (i in seq_len(nrow(ol_indx_M))) {
      temp_i <- ol_indx_M[i, 1]
      new_gene_names_list[[temp_i]] <- c(new_gene_names_list[[temp_i]],
                                         new_gene_names_temp[ol_indx_M[i, 2]])
    }
    rm(ol_indx_M, temp_i, new_gene_names_temp)
    new_gene_names_list <- lapply(new_gene_names_list, sort)
    new_gene_names <- vapply(new_gene_names_list, function(x) paste(x,
                                                                    collapse = ","), character(1))
    names(exbg) <- new_gene_names
    rm(new_gene_names, new_gene_names_list)
    rd_exons <- reduce(unlist(exbg), min.gapwidth = 0L)
    fol <- findOverlaps(rd_exons, exbg)
    split_indx <- rep(NA, length(rd_exons))
    split_indx[queryHits(fol)] <- names(exbg)[subjectHits(fol)]
    unique_exons_gene <- split(rd_exons, split_indx)
    return(unique_exons_gene)
  }
}

## A function to generate sliding window on mature RNA transcript
exonicBins <- function(exByGene,
                       binWidth = 25,
                       stepWidth = 25) {
#require(GenomicFeatures)
#Partition exons by genes
  stopifnot(stepWidth <= binWidth)

  #exByGene  <- exonsBy(txdb, by = "gene")

  tx_widths <- sum(width(exByGene))

  #Try to define the bins start always from the five prime ends of any transcripts / genes.

  bin_nums_on_tx <-
    ceiling(pmax((tx_widths - binWidth) / stepWidth, 1)) + 1 #About 7 million exome bins on hg19.

  strands_tx <- as.vector(strand(unlist(range(exByGene))))

  indx_plus <- strands_tx == "+"

  indx_minus <- strands_tx == "-"

  indx_unknown <- strands_tx == "*"

  strands_bins <- rep(strands_tx, bin_nums_on_tx)

  indx_bin_plus <- strands_bins == "+"

  indx_bin_minus <- strands_bins == "-"

  indx_bin_unknown <- strands_bins == "*"

  seqnames_bins <- rep(names(tx_widths), bin_nums_on_tx)

  bin_starts_on_tx <- vector("integer", length = sum(bin_nums_on_tx))

  bin_starts_on_tx[indx_bin_plus] <-
    unlist(lapply(bin_nums_on_tx[indx_plus], function(x)
      seq(1, stepWidth * x, by = stepWidth)), use.names = FALSE)

  bin_starts_on_tx[indx_bin_minus] <-
    unlist(mapply(
      function(x, y)
        seq(y, y - stepWidth * (x - 1), by = -1 * stepWidth),
      bin_nums_on_tx[indx_minus],
      tx_widths[indx_minus]
    ),
    use.names = FALSE) - binWidth + 1

  bin_starts_on_tx[indx_bin_unknown] <-
    unlist(lapply(bin_nums_on_tx[indx_unknown], function(x)
      seq(1, stepWidth * x, by = stepWidth)), use.names = FALSE)

  rm(bin_nums_on_tx,
     strands_tx,
     indx_plus,
     indx_minus,
     indx_unknown,
     indx_bin_plus,
     indx_bin_minus,
     indx_bin_unknown)

  bins_on_tx <- GRanges(
    seqnames = seqnames_bins,
    ranges = IRanges(start = bin_starts_on_tx,
                     width = binWidth),
    strand = strands_bins
  )

  #Trim over-hanging ends
  tx_widths <- sum(width(exByGene))

  suppressWarnings(seqlengths(bins_on_tx) <-
                     tx_widths[names(seqlengths(bins_on_tx))])

  bins_on_tx <- trim(bins_on_tx)

  bins_on_tx <- bins_on_tx[width(bins_on_tx) >= 10]

  bins_on_genome <-
    suppressWarnings(mapFromTranscripts(bins_on_tx, exByGene))

  names(bins_on_genome) <- seq_along(bins_on_genome)

  rm(bins_on_tx)

  #Removal of introns is time consuming ~ 1min.
  bins_on_genome <-
    removeIntrons(bins_on_genome, exByGene)
  
  indx <- which(sum(width(bins_on_genome)) > binWidth)
  
  if(length(indx) > 0) bins_on_genome <- bins_on_genome[-1*indx]

  return(bins_on_genome)
}

## A function to remove intronic regions from GRangesLists
removeIntrons <- function(grl,
                          exByGene){
    #Calculate intronic regions
    Introns_iranges <- gaps(ranges(exByGene))
    unlist_ebg <- unlist(exByGene)

    seq_lev <- tapply(as.vector( seqnames(unlist_ebg) ), names(unlist_ebg), function(x) x[1] )
    strand_lev <- tapply(as.vector( strand(unlist_ebg) ), names(unlist_ebg), function(x) x[1] )

    #Find the mapping between introns and bins, for only those bins that "contain" introns.
    introns_granges <- GRanges(
      seqnames = rep(seq_lev, elementNROWS(Introns_iranges)),
      ranges = unlist(Introns_iranges),
      strand = rep(strand_lev, elementNROWS(Introns_iranges))
    )

    fol <- findOverlaps(introns_granges,
                        grl,
                        type = "within")

    #Remove all the hits that are inter-genes.
    #indx_keep <- names(introns_granges)[queryHits(fol)] == gsub("\\.[0-9]*$","",names(exByGene))[grl$transcriptsHits[subjectHits(fol)]]
    indx_keep <- names(introns_granges)[queryHits(fol)] == names(exByGene)[grl$transcriptsHits[subjectHits(fol)]]
    fol <- fol[indx_keep,]

    #Split, and re-define the start and ends of those hitted bins.
    indx_Hitted_bins <-  subjectHits(fol)

    bins_contain_introns <- grl[indx_Hitted_bins]
    mcols(bins_contain_introns) <- NULL
    names(bins_contain_introns) <- indx_Hitted_bins

    #For each element within this GRanges, there is going to be one intron / multiple introns.

    introns_each_bins <- introns_granges[queryHits(fol)]
    names(introns_each_bins) <- indx_Hitted_bins

    bins_contain_introns <- c(bins_contain_introns,introns_each_bins)
    bins_contain_introns <- split(bins_contain_introns,names(bins_contain_introns))

    if(length(bins_contain_introns) == 0) {

      bins_intron_removed <- grl
      return(bins_intron_removed)

    }else{
      chunk_num = 1e5
      index_start = 1
      for(i in seq_len(ceiling( length(bins_contain_introns)/chunk_num ))) {
        Indx <- index_start:min(i*chunk_num, length(bins_contain_introns))
        bins_contain_introns[Indx] <- disjoin(bins_contain_introns[Indx])
        index_start = i*chunk_num + 1
      }

      #Remove the introns from GRanges list.
      bins_contain_introns <- unlist(bins_contain_introns)
      bins_contain_introns <- subsetByOverlaps(bins_contain_introns,
                                               introns_granges,
                                               type = "equal",invert = TRUE)
      indx_non_introns <- which( !seq_along(grl) %in% indx_Hitted_bins )

      bins_without_granges <- grl[indx_non_introns]
      mcols(bins_without_granges) <- NULL
      names(bins_without_granges) <- indx_non_introns

      bins_intron_removed <- c(bins_without_granges,bins_contain_introns)

      rm(bins_without_granges)
      rm(bins_contain_introns)

      bins_intron_removed <- bins_intron_removed[order(as.numeric(names(bins_intron_removed)))]
      names(bins_intron_removed) <- names(grl)[as.integer( names(bins_intron_removed) )]

      bins_intron_removed <- split(bins_intron_removed, names(bins_intron_removed))
      bins_intron_removed <- bins_intron_removed[order(as.numeric(names(bins_intron_removed)))]

      return(bins_intron_removed)
    }
}

## A function to flank GRangesList on the coordinate of mature RNA transcript
exonicFlank <- function(grl,
                        exByGene,
                        flankLength = 100){
  bd_on_tx <- mapToTranscripts(unlist(grl), exByGene)
  #remove names of the inner Granges (don't contain . in its name.)
  names(bd_on_tx) <- gsub("\\..*$","",names(bd_on_tx))
  bd_on_tx <- unlist( range( split(bd_on_tx, names(bd_on_tx)) ) )
  bins_on_tx <- bd_on_tx + flankLength
  rm(bd_on_tx)

  #Trim over-hanging ends
  tx_widths <- sum( width(exByGene) )
  suppressWarnings( seqlengths(bins_on_tx) <- tx_widths[names(seqlengths(bins_on_tx))] )
  bins_on_tx <- trim(bins_on_tx)
  bins_on_genome <- suppressWarnings(  mapFromTranscripts(bins_on_tx, exByGene) )
  rm(bins_on_tx)
  bins_on_genome <- trim( removeIntrons(bins_on_genome, exByGene) )
  if(is(bins_on_genome, "GRanges")) bins_on_genome <- split(bins_on_genome, names(bins_on_genome))
  return(bins_on_genome)
}

## A function to reduce GRangesList on the coordinate of mature RNA transcript
reducePeaks <- function(grl,
                        exByGene) {
  reduced_peaks_on_genome <- mapFromTranscripts( reduce( mapToTranscripts( unlist(grl) , exByGene) ), exByGene )
  names(reduced_peaks_on_genome) <- reduced_peaks_on_genome$xHits
  reduced_peaks_on_genome <- removeIntrons( reduced_peaks_on_genome, exByGene )
  if(is(reduced_peaks_on_genome, "GRanges")){
    mcols(reduced_peaks_on_genome) <- NULL
    reduced_peaks_on_genome <- split(reduced_peaks_on_genome, seq_along(reduced_peaks_on_genome))
  }
  return(reduced_peaks_on_genome)
}

## A function to calculate topology of ranges on transcripts
topologyOnTranscripts <- function(x,
                                  txdb,
                                  region_weights = c(1/3,1/3,1/3),
                                  ambiguityMethod = c("mean", "sum", "min", "max"),
                                  ignore.strand=FALSE){

  u5bytx <- fiveUTRsByTranscript(txdb)
  topology <- extractRegionRelativePosition(x,
                                            u5bytx,
                                            ambiguityMethod=ambiguityMethod,
                                            nomapValue="NA",
                                            ignore.strand=ignore.strand)*region_weights[1]
  rm(u5bytx)
  cdsbytx <- cdsBy(txdb, by = "tx")
  cdsrps <- extractRegionRelativePosition(x,
                                          cdsbytx,
                                          ambiguityMethod=ambiguityMethod,
                                          nomapValue="NA",
                                          ignore.strand=ignore.strand)
  rm(cdsbytx)
  indx <- !is.na(cdsrps)
  topology[indx] <- cdsrps[indx]*region_weights[2] + region_weights[1]
  rm(indx,cdsrps)

  u3bytx <- threeUTRsByTranscript(txdb)
  u3rps <- extractRegionRelativePosition(x,
                                         u3bytx,
                                         ambiguityMethod=ambiguityMethod,
                                         nomapValue="NA",
                                         ignore.strand=ignore.strand)
  rm(u3bytx)
  indx <- !is.na(u3rps)
  topology[indx] <- u3rps[indx]*region_weights[3] + region_weights[2] + region_weights[1]
  rm(indx,u3rps)

  return(topology)
}

## A function to calculate the relative position of ranges on transcript regions
extractRegionRelativePosition <-
function(x,
         region=NULL,
         ambiguityMethod=c("mean", "sum", "min", "max"),
         nomapValue=c("NA","0"),
         ignore.strand=FALSE){
  #require(GenomicFeatures)
  if(is(x,"GRangesList")) x <- unlist(range(x))
  x <- resize(x, width = 1, fix = "center")
  ambiguityMethod <- match.arg(ambiguityMethod)
  nomapValue <- match.arg(nomapValue)
  nomapValue <- eval(parse(text = nomapValue))

  if(is.null(region)){
    rrp_property <- rep(nomapValue, length(x))
  }else if(is(region, "GRanges")|is(region, "GRangesList")){
    if(is(region, "GRanges")){
      region_grl <- split(region, seq_along(region))
    }else{
      region <- region[elementNROWS(region) != 0]
      more_strand_region <- which(elementNROWS(runValue(strand(region))) == 1)
      region <- grl_resolve_multi_strand(region)
      names(region) <- seq_along(region)
      region_grl <- region
    }
    rrp_property <- rep(nomapValue, length(x))
    map2tx <- mapToTranscripts(x, region_grl, ignore.strand=ignore.strand)
    relpos <- start(map2tx)/sum(width(region_grl))[map2tx$transcriptsHits]
    weighted_relpos <- tapply(relpos, map2tx$xHits, eval(parse(text = ambiguityMethod[1])))
    rm(relpos)
    rrp_property[as.numeric(names(weighted_relpos))] <- weighted_relpos
    rm(map2tx, weighted_relpos)
  }else{
    stop("`region` should be either `GRanges` or `GRangesList`")
  }
  return(rrp_property)
}

## A function to resolve GRangesList mapped to multiple strands
grl_resolve_multi_strand <- function(grl){
  indx_multi <- elementNROWS(runValue(strand(grl))) > 1
  if(any(indx_multi)){
    gr_ambiguous <- unlist(grl[indx_multi])
    names(gr_ambiguous) <- paste0(names(gr_ambiguous), strand(gr_ambiguous))
    gr_resolved <- split(gr_ambiguous, names(gr_ambiguous))
    rm(gr_ambiguous)
    return(c(grl[!indx_multi],gr_resolved))
  }else{
    return(grl)
  }
}

## A helper function to reverse the read strands for 1st stranded library
readsReverse <- function(reads,
                          ...) {
  if (is(reads, "GAlignmentsList")) {
    indx_lst <- rep(seq_along(reads), elementNROWS(reads))
    reads <- unlist(reads)
    reads <- as(reads, "GRanges")
    levels( strand(reads) ) = c("-","+","*")
    reads <- split(reads, indx_lst)
    return(reads)
  } else {
    reads <- as(reads, "GRanges")
    levels( strand(reads) ) = c("-","+","*")
    return(reads)
  }
}

#'@title Silencing unwanted function output.
#'@param x any R expression.
#'
#'@return none.
#'
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

#' @title Transformation of Txdb object
#'
#' @description  Converting TxDb object into full transcript or whole genome modes.
#'
#' @param txdb a TxDb object of standard transcript annotation.
#' @param type the type of output, should be one of c("full_tx","whole_genome").
#' @return a TxDb object.
#'
#' @import GenomicFeatures
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @importFrom rtracklayer asGFF
#'
#' @name convertTxDb
#' @keywords internal
#'
convertTxDb <- function(txdb,type = c("full_transcript", "whole_genome")){
  txdb_gff <- asGFF(txdb)
  #Remove the old exon and CDS region types
  txdb_gff <- txdb_gff[txdb_gff$type != "exon"]
  txdb_gff <- txdb_gff[txdb_gff$type != "CDS"]

  if(type == "full_transcript") {
    gr_tx <- transcripts( txdb )
    mcols(gr_tx) <- DataFrame(Parent = paste0("TxID:", gr_tx$tx_id),
                              ID = NA,
                              Name = NA,
                              type = "exon")
    txdb_gff <- c(txdb_gff,gr_tx)
    rm(gr_tx)
  } else {
    if(type == "whole_genome"){
      chromosome_lengths <- seqlengths(txdb_gff)
      if(anyNA(chromosome_lengths)){
        chromosome_lengths <- tapply(end(txdb_gff),
                                     as.character(seqnames(txdb_gff)),
                                     max)
      }
      N <- length(chromosome_lengths)
      chrom_names <- names(chromosome_lengths)
      gr_plus <- GRanges(seqnames = chrom_names,
                         ranges = IRanges(start = rep(1, N ),
                                          width = chromosome_lengths),
                         strand = "+")
      gr_minus <- GRanges(seqnames = chrom_names,
                          ranges = IRanges(start = rep(1, N ),
                                           width = chromosome_lengths),
                          strand = "-")
      txdb_gff <- c(gr_plus,gr_minus,
                    gr_plus,gr_minus,
                    gr_plus,gr_minus)
      rm(gr_plus,chromosome_lengths)
      gnames <- c(paste0( chrom_names, "+" ),
                  paste0( chrom_names, "-" ))
      mcols(txdb_gff) <- DataFrame(
        type = rep(c("gene","mRNA","exon"), each = 2*N),
        ID = c(paste0("GeneID:",seq_len(2*N)),
               paste0("TxID:",seq_len(2*N)),
               rep(NA,2*N)),
        Name = c( gnames , gnames ,rep(NA,2*N) ),
        Parent = c(rep(NA,2*N),
                   paste0("GeneID:",seq_len(2*N)),
                   paste0("TxID:",seq_len(2*N)))
      )
    }
  }
  return(makeTxDbFromGRanges(txdb_gff))
}

## A function to retrieve peak annotation
makePeakAnnot <- function(peak, se, res, exbg, identity = FALSE){
  fol <- findOverlaps(peak, se, type = ifelse(identity, "equal", "any"))

  annotCols <- data.frame(pvalue = tapply(res$pvalue[subjectHits(fol)], queryHits(fol), mean, na.rm = TRUE),
                          fdr = tapply(res$padj[subjectHits(fol)], queryHits(fol), mean, na.rm = TRUE),
                          log2FC = tapply(res$log2FoldChange[subjectHits(fol)], queryHits(fol), mean, na.rm = TRUE))

  if(is.null(se$Perturbation)){
    annotCols$RPM.IP <- rpm(se, se$IP_input == "IP", fol)
    annotCols$RPM.input <- rpm(se, se$IP_input != "IP", fol)
  }else{
    annotCols$RPM.IP.Control <- rpm(se, se$IP_input == "IP" & se$Perturbation == "C", fol)
    annotCols$RPM.input.Control <- rpm(se, se$IP_input != "IP"& se$Perturbation == "C", fol)
    annotCols$RPM.IP.Treated <- rpm(se, se$IP_input == "IP" & se$Perturbation != "C", fol)
    annotCols$RPM.input.Treated <- rpm(se, se$IP_input != "IP" & se$Perturbation != "C", fol)
  }

  fol2 <- findOverlaps(peak, exbg)
  annotCols$geneID <- tapply(names(exbg)[subjectHits(fol2)], queryHits(fol2), paste, collapse = ",")

  if(is.null(se$Perturbation)){
    annotCols <- annotCols[,c("RPM.IP","RPM.input","geneID","log2FC","pvalue","fdr")]
  }else{
    colnames(annotCols)[which(colnames(annotCols)=="log2FC")] <- "diff.log2FC"
    annotCols <- annotCols[,c("RPM.IP.Control","RPM.input.Control","RPM.IP.Treated","RPM.input.Treated","geneID","diff.log2FC","pvalue","fdr")]
  }
  return(annotCols)
}

## A function to calculate RPM
rpm <- function(se, col_indx, fol){
  count <- rowSums(cbind(assay(se)[, col_indx]))
  RPM <- count/sum(count)*1e6
  RPM <- tapply(RPM[subjectHits(fol)], queryHits(fol), sum)
  RPM <- round(RPM, 2)
  return(RPM)
}

## A function to save peak statistics
savePeak <- function(peak, save_dir = "exomePeak2_output", file_name = c("peaks","diffPeaks")){
  if (length(peak) == 0){
  message("No significant peaks detected, result unsaved.")
  }else{
  if (!dir.exists(save_dir)) dir.create(save_dir)
  mcols(peak)$score <- -log10(mcols(peak)$pvalue)
  mcols(peak)$score[is.na(mcols(peak)$score)] <- 0
  export(peak, file.path(save_dir, paste0(file_name, ".bed")))
  Tbl <- read.table(file.path(save_dir, paste0(file_name, ".bed")), header = FALSE, sep = "\t")
  Tbl <- Tbl[,-1*c(5,7,8,9)]
  colnames(Tbl) <- c("chr","chromStart",
                     "chromEnd","name",
                     "strand","blockCount",
                     "blockSizes", "blockStarts")
  Tbl <- cbind(Tbl, mcols(peak))
  rownames(Tbl) <- NULL
  write.csv(Tbl, file.path(save_dir, paste0(file_name, ".csv")))
  saveRDS(peak, file.path(save_dir, paste0(file_name, ".rds")))
  }
}

## Motif sampler
sampleSequence <- function(motif, region, sequence, fixed = FALSE, replace = FALSE){
  #require(BSgenome)
  #require(GenomicFeatures)
  stopifnot(is(region, "GRangesList")|is(region, "GRanges"))
  if(is(region, "GRangesList")) region <- unlist(region)
  region <- reduce(region)

  region_dnass <- getSeq(x=sequence,
                         names=seqnames(region),
                         start=start(region),
                         end=end(region),
                         strand=strand(region),
                         as.character=FALSE)

  indx <- paste0("reg_", seq_along(region))
  regions_GRL <- split(region, indx)
  regions_GRL <- regions_GRL[indx]
  rm(indx)
  vmp <- vmatchPattern(motif, region_dnass, fixed = fixed)
  rm(region_dnass)
  vmp_gr <- GRanges(seqnames = rep(names(regions_GRL), elementNROWS(vmp)), ranges = unlist(vmp))
  rm(vmp)
  motif_on_regions <- mapFromTranscripts(vmp_gr,regions_GRL)
  rm(vmp_gr, regions_GRL)
  mcols(motif_on_regions) = NULL
  seqlengths(motif_on_regions) <- seqlengths(region)
  return(motif_on_regions)
}



