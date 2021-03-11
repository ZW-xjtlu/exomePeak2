#'@title a function to sort the summarizedExomePeak by its order on genome
#'@param sep a summarizedExomePeak object.
#'
#'@return the summarizedExomePeak with a given order.
#'@keywords internal
#'@importFrom BiocGenerics order
#'
sort_sep<-function(sep){
  order_indx <- order(unlist(range(rowRanges(sep))))
  sep <- sep[order_indx,]
  exomePeak2Results(sep) <- exomePeak2Results(sep)[order_indx,]
  return(sep)
}
