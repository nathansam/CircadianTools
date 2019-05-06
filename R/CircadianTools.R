#' CircadianTools: Collection of Tools for Detecting Rhythmic Genes
#'
#'Provides functions to detect rhythmicity using Cosinor models or turning point analysis.
#' @section General Tools
#' \code{\link{basicplot}}: Plots activity data as points and mean data as lines
#' @section Cosinor Functions:
#' \code{\link{cosinoranalysis}}: Fits cosinor models to transcriptomics data and plots the best-fitting models using ggplot2. \cr \cr
#' \code{\link{cosinorplot}}: Fits a cosinor model to a given gene in a given dataset and plots the model. \cr \cr
#' @section Turning Point Analysis:
#' \code{\link{turningplot}}: Fits a spline to a given gene in a given dataset. Finds the turning points. Plots the turning points and spline.
#' @section Utility Functions:
#' \code{\link{ggplot.cosinor.lm}}: Adapted from the Cosinor package by Michael Sachs. Given a cosinor.lm model fit, generate a plot of the data with the fitted values.
#' \code{\link{maketimevector}}: Produces a vector of time values for the gene activity readings. \cr \cr
#' \code{\link{geneclean}}: Removes columns and rows which show no gene activity over time. \cr \cr
#'
#' @docType package
#' @name CircadianTools
NULL
