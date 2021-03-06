% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CircadianTools.R
\docType{package}
\name{CircadianTools}
\alias{CircadianTools}
\title{CircadianTools: A Collection of Tools for Detecting Rhythmic Genes}
\description{
Allows analysis of rhythmic genes to be easily carried out on transcriptomics
data using R. Designed to be as flexible as possible such as by allowing an
unequal number of replicates. \cr \cr
------------------------------------------------------------------------------
}
\section{Simple Plots}{

\code{\link{BasicPlot}}: Plots activity data as points and average activity
 as lines \cr \cr

\code{\link{CompPlot}}: Plots two genes from a gene activity dataset. \cr \cr

\code{\link{DatasetPlot}}: Saves plots of all genes in a dataset. \cr \cr

\code{\link{TurningPlot}}: Fits a spline to a given gene in a given dataset.
 Finds the turning points and then plots the turning points and spline.
 \cr \cr
------------------------------------------------------------------------------
}

\section{Clustering}{


\code{\link{AgglomClustering}}: Applies hierarchical clustering, clustering to a
 transcriptomics dataset and appends a cluster column to this dataset for all
 genes. \cr \cr

\code{\link{AgglomParamSelection}}: Runs hierarchical clustering with
 differing numbers of partitions and returns validation metrics. \cr \cr

\code{\link{ClusterCenterGenerator}}: Finds the center of every cluster in a
 dataset \cr \cr

\code{\link{ClusterCorDatasetPlot}}: Uses \link{ClusterCorPlot} to plot all
 of the clusters generated by a clustering method when absolute
 Pearson's correlation was used as a distance measure. \cr \cr

 \code{\link{ClusterCorPlot}}: Plots the activity level for a cluster
 generated by using absolute Pearson's correlation as a distance measure.
 Plots positively and negatively correlated genes as two different lines.
 \cr \cr

\code{\link{ClusterDatasetPlot}}: Plots the mean and error bars for all
 clusters across time. \cr \cr

\code{\link{ClusterParamSelection}}: Calculates validation metrics for
 different clustering methods and different numbers of partitions.
 The validation metrics are plotted. \cr \cr

\code{\link{ClusterPlot}}: Plots the mean and error bars for the genes in a
 cluster across time. \cr \cr

\code{\link{ClusterSpread}}: Shows how many genes are in each cluster after
 clustering has been applied. \cr \cr

\code{\link{ClusterText}}: Takes a dataframe of clusters and stores the name
 of all genes in a text file. The row number deontes the cluster number. \cr
 \cr

\code{\link{ClusterTimeProfile}}: Provides a dataframe of median values at
 each time point for each cluster. \cr \cr

\code{\link{CommonSingletonFinder}}: Finds the genes which belong to common
 singleton clusters in two different clustered datasets. \cr \cr

 \code{\link{DendogramDatasetPlot}}: Plots the dendogram for every cluster in
  a clustered dataset. \cr \cr

\code{\link{DendogramPlot}}: Plots the dendogram for a cluster in a
 clustered dataset \cr \cr

\code{\link{DianaClustering}}: Applies DIANA (DIvisive ANAlysis) clustering
 to a transcriptomics dataset and appends a cluster column to this dataset
 for all genes. \cr \cr

\code{\link{DianaParamSelection}}:  Runs DIANA (DIvisive ANAlysis) clustering
 with differing numbers of partitions and returns validation metrics. \cr \cr

\code{\link{FindClusterDistanceQuantiles}}: Finds The distances between the
 center of each cluster and the centers of all the other clusters. \cr \cr

\code{\link{FindClusterMedian}}: Finds the center of a cluster \cr \cr

\code{\link{FindClusterQuantile}}: Finds The distances between the center of
 a cluster and the centers of all other clusters. \cr \cr

 \code{\link{MDSPlot}}: Applies multidimensional scaling to a clustered
  transcriptomics dataset to reduce the clusters to two dimensions and
   then plots the clusters.  \cr \cr

\code{\link{PamClustering}}: Applies PAM (Partitioning around Medoids)
 clustering to a transcriptomics dataset and appends a cluster column to this
 dataset for all genes. \cr \cr

\code{\link{PamParamSelection}}: Runs PAM with differing numbers of
 partitions and returns validation metrics. \cr \cr

\code{\link{QuantilePlots}}: Finds the quartiles for intercluster distances
 and plots these distances as a set of histograms \cr \cr

\code{\link{SingletonNameFinder}}: Finds the genes which belong to
 singleton clusters. \cr \cr
 ----------------------------------------------------------------------------
}

\section{Correlation}{

\code{\link{CorAnalysis}}: Ranks correlation between a given gene and all
 other genes in a dataset. Plots both the given gene and highly correlated
 genes for a given correlation value. \cr \cr

\code{\link{CorAnalysisCluster}}: Correlates the average activity of a
 cluster with the average activity of every other cluster. \cr \cr

\code{\link{CorAnalysisClusterDataset}}: Correlates the average activity of
 each cluster with every other cluster in a dataset. \cr \cr

\code{\link{CorAnalysisDataset}}: Correlates every gene in a dataset with
 every other gene in the same dataset. Allows a timelag between genes to be
 correlated. \cr \cr

\code{\link{CorAnalysisPar}}: Parallel Implementation of \code{CorAnalysis}.
 \cr \cr


\code{\link{CorSignificantPlot}}: Prints or saves the genes found to be
 most significant by \code{CorAnalysis}. \cr \cr
 ----------------------------------------------------------------------------
}

\section{Cosinor}{

\code{\link{CosinorAnalysis}}: Fits cosinor models to transcriptomics data
 and plots the best-fitting models using ggplot2. \cr \cr

\code{\link{CosinorAnalysisPar}}: Parallel Implementation of
\code{CosinorAnalysis}. \cr \cr

\code{\link{CosinorPlot}}: Fits a cosinor model to a given gene in a given
 dataset and plots the model. \cr \cr

\code{\link{CosinorResidualDatasetPlot}}: Fits cosinor models and plots the
 residuals for multiple genes in a dataset \cr \cr

\code{\link{CosinorResidualPlot}}: Fits a cosinor model to a gene and plots
 the residuals \cr \cr

\code{\link{CosinorSignificantPlot}}: Prints or saves the genes found to be
 most significant by \code{cosinoranalysis}. \cr \cr

\code{\link{MultiCosinorTest}}: Fits a cosinor model and carries out ANOVA
 using raw coefficients. Then fits a cosinor model with additonal sine and
 cosine terms with a different period. ANOVA tests are carried out on the
 more complex model as well as directly comparing the two models. \cr \cr
 ----------------------------------------------------------------------------
}

\section{Cytoscape}{

\code{\link{CytoscapeFile}}: Converts a correlation dataframe object into a
 format suitable for cytoscape and saves as a csv file. \cr \cr

\code{\link{CytoscapeFilter}}: Reduces the size of a file intended for
Cytoscape by filtering out the genes/clusters which are not correlated
\cr \cr
-----------------------------------------------------------------------------
}

\section{Fasta Files}{

\code{\link{ContigGen}}: Finds all unique contig IDs in a transcriptomics
 dataset \cr \cr

\code{\link{FastaSub}}: Creates a fasta file from only certain sequences in
 another fasta file \cr \cr
 ----------------------------------------------------------------------------
}

\section{Filtering}{

\code{\link{AnovaFilter}}: Filters a gene activity dataframe via ANOVA.
\cr \cr

\code{\link{CombiFilter}}: Filters a transcriptomics dataset by using
 \code{Zerofilter}, \code{Anovafilter} and \code{Sizefilter}. \cr \cr

\code{\link{SizeFilter}}: Filters the genes with the smallest range from a
 transcriptomics dataset. \cr \cr

\code{\link{TFilter}}: Applies a filter where a t-test is carried out on gene
 activity levels between time points. \cr \cr

\code{\link{ZeroFilter}}: Filters a transcriptomics dataset such that there
 is a minimum number of non-zero activity readings for each gene in the
 resultant dataset. \cr \cr
 ----------------------------------------------------------------------------
}

\section{RAIN}{

\code{\link{RainAnalysis}}: Carries out RAIN analysis on a transcriptomics
 dataset. \cr \cr

\code{\link{RainSignificantPlot}}: Prints or saves plots of the genes found
 to be most significant by \code{rainanalysis}. \cr \cr
-----------------------------------------------------------------------------
}

\section{Utility Functions}{

\code{\link{AbsCorDist}}: Calculates a distance matrix based on the distance
 measure of: \deqn{1 - |cor(x, y)|} \cr \cr

\code{\link{ActivitySelect}}: Returns gene activity by either gene name or
 row number \cr \cr

\code{\link{FileConflict}}: Checks if a file which will be created already
 exists and, if necessary asks the user if this file should be overwritten.
 \cr \cr

\code{\link{GeneClean}}: Removes columns and rows which show no gene activity
 over time. \cr \cr

\code{\link{GeneRange}}: Finds the range of gene activity for each gene in a
 dataframe. The median for the replicates is used for each time point.
 \cr \cr

\code{\link{GeneScale}}: Centers/scales every gene in a transcriptomics
 dataset. \cr \cr

\code{\link{GeneSub}}: Takes an object where the first column is genenames
 (IE a column of known Circadian genes) and subsets from a dataset containing
 activity readings for these genes. \cr \cr

\code{\link{ggplot.cosinor.lm}}: Adapted from the Cosinor package by
Michael Sachs. Given a cosinor.lm model fit, generate a plot of the data with
the fitted values.\cr \cr

\code{\link{MakeTimevector}}: Produces a vector of time values for the gene
 activity readings. \cr \cr

\code{\link{MedList}}: Provides a dataframe of median values at each
 time point for each gene from a transcriptomics dataset. \cr \cr

\code{\link{TAnalysis}}: A t.test is carried out on gene
 activity levels between time points and the number of significant increases
 & decreases is returned. \cr \cr
}

