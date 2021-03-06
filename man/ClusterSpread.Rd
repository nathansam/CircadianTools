% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_functions.R
\name{ClusterSpread}
\alias{ClusterSpread}
\title{ClusterSpread:}
\usage{
ClusterSpread(cluster.dataset)
}
\arguments{
\item{cluster.dataset}{A transcriptomics dataset where the final column
details the cluster the gene belongs to. First column should be gene names.
All remaining columns should be expression levels.}
}
\value{
A dataframe object. The first column is the cluster number. Second
 column is how many genes belong to that cluster.
}
\description{
Shows how many genes are in each cluster after clustering has
 been applied.
}
\examples{
pam.df <- PamClustering(Laurasmappings, k = 10)
clusterstats<-ClusterSpread(pam.df)
}
