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
