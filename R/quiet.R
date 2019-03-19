#'@title Silencing unwanted function output.
#'@param x any R expression.
#'
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
