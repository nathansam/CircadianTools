% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simple_plot_functions.R
\name{CompPlot}
\alias{CompPlot}
\title{CompPlot:}
\usage{
CompPlot(gene1, gene2, dataset, save = FALSE, points = TRUE)
}
\arguments{
\item{gene1}{The name of the first gene to be plotted. Must be a string.}

\item{gene2}{The name of the second gene to be plotted. Must be a string.}

\item{dataset}{A transcriptomics dataset. First columns should be gene names.
All other columns should be expression levels.}

\item{save}{Logical. If TRUE, saves plots to working directory.
Defaults to FALSE.}

\item{points}{Logical. If FALSE then each observation will not be plotted
as points.}
}
\value{
Returns or saves a ggplot2 object.
}
\description{
Plots two genes from a gene activity dataset
}
\examples{
CompPlot('comp100000_c0_seq3', 'comp100002_c0_seq2', Laurasmappings)

}
