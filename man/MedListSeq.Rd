% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{MedListSeq}
\alias{MedListSeq}
\title{MedListSeq}
\usage{
MedListSeq(dataset)
}
\arguments{
\item{dataset}{A transcriptomics dataset. First columns should be gene names.
All other columns should be expression levels.}
}
\description{
Provides a dataframe of median values at each time point for
 each gene from a transcriptomics dataset.
}
\examples{
mediandf <- MedListSeq(Laurasmappings)

}
