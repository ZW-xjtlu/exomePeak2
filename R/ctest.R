#'@title Exact poisson test on the ratio between IP and input counts.
#'@param IP_count a numeric vector for the counts of IP sample.
#'@param input_count a numeric vector for the counts of input sample.
#'@param IP_sizeFactor a numeric vector for the size factors of IP sample.
#'@param input_sizeFactor a numeric vector for the size factors of input sample.
#'@param fold the fold change under the null hypothesis; default = 1.
#'@param alpha the alpha level used for determining the optimal independent filtering in \code{filtered_R}; default = 0.05.
#'@details C-tests will be conducted between each entries of the IP and input vector.
#'
#' The following statistical model is assumed on the data:
#'
#' IP ~ Poisson(mean_IP*IP_sizeFactor)
#'
#' input ~ Poisson(mean_input*input_sizeFactor)
#'
#' The one-sided statistical test (C-test) is conducted with the following hypothesis pair:
#'
#' \deqn{H_0: mean_IP/mean_input <= fold}
#'
#' \deqn{H_1: mean_IP/mean_input > fold}
#'
#' The exact p-values will be generated using binomial test, check \code{\link{poisson.test}} for more details.
#'
#' The p-values are adjusted with "BH" method with an independent filtering slected by function \code{filtered_R} in package \code{genefilter};
#'
#' the filter statistics used is the normalized read abundence.
#'
#'@import SummarizedExperiment
#'@import genefilter
#'
#'@return A list include adjusted p values (with method fdr) and the corresponding log2 fold changes.

ctest <-  function(IP_count,
                   input_count,
                   IP_sizeFactor,
                   input_sizeFactor,
                   fold=1,
                   alpha = 0.05) {

  # input check
  if (length(IP_count) != length(input_count)) { stop("The IP and INPUT of ctest must be of the same length.", call. = TRUE, domain = NULL) }

  # replace 0 with 1
  IP_count = pmax(IP_count,1)
  input_count = pmax(input_count,1)

  # calculate p
  a=IP_sizeFactor*fold
  b=input_sizeFactor
  p=a/(a+b)

  # get total observation
  total=IP_count+input_count

  # p value and adjusted p values with independent filtering
  pvalue = pbinom(IP_count-1, total, p, lower.tail = FALSE, log.p = FALSE)

  theta = seq(from=0, to=0.8, by=0.02)

  rejBH = filtered_R(alpha=alpha,
                     filter=IP_count/IP_sizeFactor + input_count/input_sizeFactor,
                     test=pvalue,
                     theta=theta,
                     method="BH")

  pBH = filtered_p(filter=IP_count/IP_sizeFactor + input_count/input_sizeFactor,
                   test=pvalue,
                   theta=theta[which.max(rejBH)],
                   method="BH")

  # fold enrichment
  log2.fc = log2((pmax(1,IP_count)/IP_sizeFactor)/(pmax(1,input_count)/input_sizeFactor))

  # output result
  PW=list(pvalue = pvalue,fdr = pBH,log2FC=log2.fc)

  return(PW)
}
