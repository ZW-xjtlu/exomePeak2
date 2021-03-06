% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/plotExonLength.R
\name{plotExonLength-methods}
\alias{plotExonLength-methods}
\alias{plotExonLength,SummarizedExomePeak-method}
\alias{plotExonLength}
\title{Method plotExonLength}
\usage{
plotExonLength(
  sep,
  txdb = NULL,
  save_pdf_prefix = NULL,
  include_control_regions = TRUE,
  save_dir = "."
)

\S4method{plotExonLength}{SummarizedExomePeak}(
  sep,
  txdb = NULL,
  save_pdf_prefix = NULL,
  include_control_regions = TRUE,
  save_dir = "."
)
}
\arguments{
\item{sep}{a \code{\link{SummarizedExomePeak}} object.}

\item{txdb}{a \code{\link{TxDb}} object containing the transcript annotation.}

\item{save_pdf_prefix}{a \code{character} if provided, a pdf file with the given name will be saved under the current working directory.}

\item{include_control_regions}{a \code{logical} for whether to include the control regions or not; Default \code{= TRUE}.}

\item{save_dir}{a \code{character} for the directory to save the plot; Default \code{= "."}.}
}
\value{
a ggplot object
}
\description{
This function plot the distribution of the exon length for peaks containing exons.
}
\details{
If the SummarizedExomePeaks object contains LFC statistics, the significantly modified peaks
with IP to input log2FC > 0 and GLM Wald test padj < .05 will be plotted .

If the SummarizedExomePeaks object contains interactive LFC statistics, both the hyper modification
and hypo modification peaks with GLM Wald test p values < .05 will be plotted.
}
\examples{

### Make TxDb object from the gff file
library(GenomicFeatures)
GENE_ANNO_GTF = system.file("extdata", "example.gtf", package="exomePeak2")

txdb <- makeTxDbFromGFF(GENE_ANNO_GTF)

### Load the example SummarizedExomPeak object
f1 = system.file("extdata", "sep_ex_mod.rds", package="exomePeak2")

sep <- readRDS(f1)

### Visualize the linear relationships between GC content and normalized reads count under different regions
plotExonLength(sep,txdb)

}
