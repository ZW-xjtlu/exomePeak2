#'@title a function to split a GRanges into GRangesList by its names.
#'@param x a list like object.
#'
#'@return the list split by its names.
#'@keywords internal
#'
split_by_name <- function(x) {
  return(split(x, names(x)))
}
