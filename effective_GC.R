#GCAPC
EGCC <- function(
                 genome = "hg19",
                 gctype = c("ladder", "tricube", "naive"))
{
  genome <- getBSgenome(genome)
  model <- match.arg(model)
  gctype <- match.arg(gctype)
  bdw <- bdwidth[1]
  halfbdw <- floor(bdw/2)
  if ((length(bdwidth) < 2 & is.null(flank)) || (!is.null(flank) &
                                                 length(bdwidth) < 1))
    stop("Parameter 'bdwidth' or 'flank' error.\n")
  if (is.null(flank)) {
    pdwh <- bdwidth[2]
    flank <- pdwh - bdw + halfbdw
  }
  else {
    pdwh <- flank + bdw - halfbdw
  }
  if (length(sampling) < 2 || sampling[1] <= 0 || sampling[1] >
      1 || sampling[2] < 1)
    stop("Parameter 'sampling' error.\n")
  sampling[2] <- floor(sampling[2])
  if (class(supervise) != "GRanges")
    stop("Parameter 'supervise' error.\n")
  if (sum(gcrange < 0) > 0 || sum(gcrange > 1) > 0 || sum(is.na(gcrange)) >
      0)
    stop("Parameter 'gcrange' error.\n")
  if (mu0 <= 0 || mu1 <= 0 || mu0 >= mu1)
    stop("Parameter 'mu0' or 'mu1' error.\n")
  if (model == "nbinom" && (theta1 <= 0 || theta0 <= 0))
    stop("Parameter 'theta0' or 'theta1' error in nbinom model.\n")
  if (p <= 0 || p >= 0.5)
    stop("'p' must be in (0,0.5).\n")
  if (converge <= 0 || converge > 0.1)
    stop("'converge' must be in (0,0.1].\n")
  if (gctype == "ladder") {
    weight <- c(seq_len(flank), rep(flank + 1, bdw), rev(seq_len(flank)))
    weight <- weight/sum(weight)
  }
  else if (gctype == "tricube") {
    w <- flank + halfbdw
    weight <- (1 - abs(seq(-w, w)/w)^3)^3
    weight <- weight/sum(weight)
  }
  cat("Starting to estimate GC effects\n")
  seqs <- sapply(coverage$fwd, length)
  seqs <- floor(seqs/bdw - 10) * bdw
  binl <- sapply(seqs, function(i) length(seq(1 + bdw * 10,
                                              i, bdw)))
  starts <- unlist(lapply(seqs, function(i) seq(1 + bdw * 10,
                                                i, bdw)))
  ends <- unlist(lapply(seqs, function(i) seq(bdw * 11, i,
                                              bdw)))
  chrs <- rep(names(seqs), times = binl)
  region <- GRanges(chrs, IRanges(start = starts, end = ends))
  rm(seqs, binl, starts, ends, chrs)
  fo <- unique(subjectHits(findOverlaps(supervise, region)))
  if (length(fo) < 1000 || length(region) - length(fo) < 1000) {
    cat("...... Too few/much ranges in 'supervise'.", "Selecting random sampling\n")
    sampidx <- sort(as.integer(sapply(seq_len(sampling[2]),
                                      function(x) sample.int(length(region), floor(length(region) *
                                                                                     sampling[1])))))
    samptype <- "random"
  }
  else {
    cat("...... Selecting supervised sampling\n")
    num1 <- floor(length(fo) * sampling[1])
    num2 <- floor((length(region) - length(fo)) * sampling[1])
    num1 <- min(num1 * ceiling(0.1/(num1/num2)), length(fo))
    sampidx <- sort(as.integer(sapply(seq_len(sampling[2]),
                                      function(x) c(sample(fo, num1), sample(setdiff(seq_along(region),
                                                                                     fo), num2)))))
    samptype <- "supervised"
  }
  region <- region[sampidx]
  cat("......... Sampling", sampling[1] * 100, "percent of regions by",
      sampling[2], "times, total", length(region), "regions\n")
  regionsp <- resize(split(region, seqnames(region)), pdwh)
  cat("...... Counting reads\n")
  rcfwd <- unlist(viewSums(Views(coverage$fwd,
                                 ranges(shift(regionsp,-flank)))))

  rcrev <- unlist(viewSums(Views(coverage$rev,
                                 ranges(shift(regionsp,halfbdw)))))
  rm(regionsp, sampidx, fo)
  cat("...... Calculating GC content with flanking",
      flank,
      "\n")
  rwidth <- width(region[1])
  nr <- shift(resize(region, rwidth + flank * 2), -flank)
  seqs <- getSeq(genome, nr)
  gcpos <- startIndex(vmatchPattern("S", seqs, fixed = "subject"))
  gc <- round(sapply(gcpos, function(x) sum(weight[x])), 3)
  rm(rwidth, nr, seqs, gcpos, region)
  return(gc)
}
