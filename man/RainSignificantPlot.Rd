% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rain_functions.R
\name{RainSignificantPlot}
\alias{RainSignificantPlot}
\title{RainSignificantPlot :}
\usage{
RainSignificantPlot(results, dataset, number = 10, print = TRUE, save = FALSE)
}
\arguments{
\item{results}{A dataframe generated by \code{RainAnalysis}}

\item{dataset}{A transcriptomics dataset which was used with
\link{RainAnalysis} to generate the results dataframe. First columns should
be gene names. All other columns should be expression levels.}

\item{number}{The number of most significant genes printed or saved}

\item{print}{Logical. If TRUE renders significant genes in the plot viewer.
Defaults to TRUE}

\item{save}{Logical. If TRUE, saves plots to working directory.
Defaults to FALSE.}
}
\value{
Prints or saves ggplot2 object(s).
}
\description{
Prints or saves plots of the genes found to be most significant
by \link{RainAnalysis}
}
\examples{
results <- RainAnalysis(Laurasmappings, period = 24)
RainSignificantPlot(results,Laurasmappings, save=TRUE, number=15)

}
