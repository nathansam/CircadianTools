% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_functions.R
\name{FindClusterMedian}
\alias{FindClusterMedian}
\title{FindClusterMedian}
\usage{
FindClusterMedian(cluster.no, cluster.dataset, nthreads = NULL)
}
\arguments{
\item{cluster.no}{The number which identifies the cluster}

\item{cluster.dataset}{A transcriptomics dataset where the final column
details the cluster the gene belongs to. First column should be gene names.
All remaining columns should be expression levels.}

\item{nthreads}{Number of processor threads for the process. If not specifed
then the maximum number of logical cores are used.}
}
\description{
Finds the center of a cluster by finding the median time value
 for each gene and then calculates the median activity for each of these
 time points across an entire cluster.
}
\examples{
pam.df <- PamClustering(Laurasmappings, k = 10, nthreads = 2)
FindClusterMedian(2, cluster.dataset = pam.df, nthreads = 2)
}
