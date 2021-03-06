% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_functions.R
\name{ClusterPlot}
\alias{ClusterPlot}
\title{ClusterPlot :}
\usage{
ClusterPlot(
  clusterno,
  cluster.dataset,
  nthreads = NULL,
  print = TRUE,
  save = FALSE,
  path = NULL
)
}
\arguments{
\item{clusterno}{The number which identifies the cluster}

\item{cluster.dataset}{A transcriptomics dataset where the final column
details the cluster the gene belongs to. First column should be gene names.
All remaining columns should be expression levels.}

\item{nthreads}{Number of processor threads for the process. If not specified
then the maximum number of logical cores are used.}

\item{print}{Logical. If TRUE renders significant genes in the plot viewer.
Defaults to TRUE}

\item{save}{Logical. If TRUE, saves plots. Defaults to FALSE.}

\item{path}{The directory to be used for saving plots to. Uses the working
directory by default. Not used if save=FALSE}
}
\value{
Prints or saves a ggplot2 object.
}
\description{
Plots the mean and error bars for the genes in a cluster across
 time
}
\examples{
pam.df <- PamClustering(Laurasmappings, k = 10, nthreads = 2)
ClusterPlot(2, pam.df, nthreads = 2)
}
