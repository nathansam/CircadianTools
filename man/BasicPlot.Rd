% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simple_plot_functions.R
\name{BasicPlot}
\alias{BasicPlot}
\title{BasicPlot:}
\usage{
BasicPlot(
  genename,
  dataset,
  timelag = 0,
  method = "median",
  points = TRUE,
  print = TRUE,
  save = FALSE,
  path = NULL
)
}
\arguments{
\item{genename}{The name of a gene intended for plotting. Must be a string.}

\item{dataset}{A transcriptomics dataset. First columns should be gene names.
All other columns should be expression levels.}

\item{timelag}{Shifts the plot to earlier in time.}

\item{method}{How should the average activity for each time point be
calculated? 'median' or 'mean'}

\item{points}{Logical. If FALSE then each observation will not be plotted
as points.}

\item{print}{Logical. If TRUE renders plot in the plot viewer. Defaults to
TRUE}

\item{save}{Logical. If TRUE, saves plot to working directory. Defaults to
FALSE.}

\item{path}{The directory to be used for saving plots to. Uses the working
directory by default. Not used if save = FALSE}
}
\value{
A ggplot2 object
}
\description{
Plots activity data as points and average activity as lines
}
\examples{
BasicPlot('comp100026_c0_seq2',Laurasmappings)

}
