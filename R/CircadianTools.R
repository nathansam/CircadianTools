#' CircadianTools: A Collection of Tools for Detecting Rhythmic Genes
#'
#' Allows analysis of rhythmic genes to be easily carried out on transcriptomics data using R. Designed to be as flexible as possible such as by allowing an unequal number of replicates. \cr \cr
#' @section General Tools:
#' \code{\link{basicplot}}: Plots activity data as points and average activity as lines \cr \cr
#' \code{\link{compplot}}: Plots two genes from a gene activity dataset. \cr \cr
#' \code{\link{datasetplot}}: Saves plots of all genes in a dataset. WARNING! Don't run on a large dataset! Intended for a filtered dataset \cr \cr
#' @section Clustering:
#' \code{\link{clusteroverview}}: Plots the mean and error bars for all clusters across time. \cr \cr
#' \code{\link{clusterplot}}: Plots the mean and error bars for the genes in a cluster across time. \cr \cr
#' \code{\link{clusterspread}}: Shows how many genes are in each cluster after clustering has been applied. \cr \cr
#' \code{\link{hclustering}}: Applies hierarchical clustering, clustering to a transcriptomics dataset and appends a cluster column to this dataset for all genes. \cr \cr
#' \code{\link{pamclustering}}: Applies PAM (Partitioning around Medoids) clustering to a transcriptomics dataset and appends a cluster column to this dataset for all genes. \cr \cr
#' @section Correlation:
#' \code{\link{coranalysis}}: Ranks correlation between a given gene and all other genes in a dataset. Plots both the given gene and highly correlated genes for a given correlation value. \cr \cr
#' \code{\link{coranalysisdataset}}: Correlates every gene in a dataset with every other gene in the same dataset. Allows a timelag between genes to be correlated. \cr \cr
#'  \code{\link{coranalysispar}}: Parallel Implementation of \code{coranalysis}. \cr \cr
#' \code{\link{corsignificantplot}}: Prints or saves the genes found to be most significant by \code{cosinoranalysis}. \cr \cr
#' @section Cosinor Functions:
#' \code{\link{cosinoranalysis}}: Fits cosinor models to transcriptomics data and plots the best-fitting models using ggplot2. \cr \cr
#' \code{\link{cosinoranalysispar}}: Parallel Implementation of \code{cosinoranalysis}. \cr \cr
#' \code{\link{cosinorplot}}: Fits a cosinor model to a given gene in a given dataset and plots the model. \cr \cr
#' \code{\link{cosinorsignificantplot}}: Prints or saves the genes found to be most significant by \code{cosinoranalysis}. \cr \cr
#' @section Filtering:
#' \code{\link{anovafilter}}: Filters a gene activity dataframe via ANOVA. \cr \cr
#' \code{\link{combifilter}}: Filters a transcriptomics dataset by using \code{zerofilter}, \code{anovafilter} and \code{sizefilter}. \cr \cr
#' \code{\link{sizefilter}}: Filters the genes with the smallest range from a transcriptomics dataset. \cr \cr
#' \code{\link{tfilter}}: Experimental! Applies a filter where a t.test is carried out on gene activity levels between time points. \cr \cr
#' \code{\link{zerofilter}}: Filters a transcriptomics dataset such that there is a minimum number of non-zero activity readings for each gene. \cr \cr
#' @section RAIN Functions:
#' \code{\link{rainanalysis}}: Carries out RAIN analysis on a gene dataset.  \cr \cr
#' \code{\link{rainsignificantplot}}: Prints or saves the genes found to be most significant by \code{rainanalysis}. \cr \cr
#' @section Turning Point Functions:
#' \code{\link{turningplot}}: Fits a spline to a given gene in a given dataset. Finds the turning points. Plots the turning points and spline.\cr \cr
#' @section Utility Functions:
#' \code{\link{geneclean}}: Removes columns and rows which show no gene activity over time. \cr \cr
#' \code{\link{generange}}: Finds the range of gene activity for each gene in a dataframe. The median for the replicates is used for each time point. \cr \cr
#' \code{\link{genesub}}: Takes an object where the first column is genenames (IE a column of known Circadian genes) and subsets from a dataset containing activity for these genes. \cr \cr
#' \code{\link{ggplot.cosinor.lm}}: Adapted from the Cosinor package by Michael Sachs. Given a cosinor.lm model fit, generate a plot of the data with the fitted values.\cr \cr
#' \code{\link{maketimevector}}: Produces a vector of time values for the gene activity readings. \cr \cr
#' \code{\link{medlist}}: Provides a dataframe of median values at each time point for each gene from a transcriptomics dataset. \cr \cr
#' \code{\link{tanalysis}}: Experimental! A t.test is carried out on gene activity levels between time points and the number of significant increases & decreases is returned. \cr \cr
#' @docType package
#' @name CircadianTools
NULL
