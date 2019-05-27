#' CircadianTools: Collection of Tools for Detecting Rhythmic Genes
#'
#' Allows Cosinor Models and Turning Point Analysis to be easily carried out on transcriptomics data using R. Designed to be as flexible as possible such as by allowing an unequal number of replicates. \cr \cr
#' @section General Tools:
#' \code{\link{basicplot}}: Plots activity data as points and mean data as lines \cr \cr
#' @section Cosinor Functions:
#' \code{\link{cosinoranalysis}}: Fits cosinor models to transcriptomics data and plots the best-fitting models using ggplot2. \cr \cr
#' \code{\link{cosinorplot}}: Fits a cosinor model to a given gene in a given dataset and plots the model. \cr \cr
#' @section Turning Point Analysis:
#' \code{\link{turningplot}}: Fits a spline to a given gene in a given dataset. Finds the turning points. Plots the turning points and spline.\cr \cr
#' @section Utility Functions:
#' \code{\link{ggplot.cosinor.lm}}: Adapted from the Cosinor package by Michael Sachs. Given a cosinor.lm model fit, generate a plot of the data with the fitted values.\cr \cr
#' \code{\link{maketimevector}}: Produces a vector of time values for the gene activity readings. \cr \cr
#' \code{\link{geneclean}}: Removes columns and rows which show no gene activity over time. \cr \cr
#'
#' @docType package
#' @name CircadianTools
NULL
