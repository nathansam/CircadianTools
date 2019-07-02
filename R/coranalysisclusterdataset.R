#' coranalysisclusterdataset:
#'
#' @description Correlates the average activity of each cluster with every other cluster in a dataset.
#' @param cluster.dataset A transcriptomics dataset where the final column details the cluster the gene belongs to. First column should be gene names. All remaining columns should be expression levels.
#' @param lag Setting any value other than 0 allows a gene to be correlated with lagged genes in the dataset. The number denotes the number of timesteps to lag by.
#' @param nthreads The number of threads to be used for parallel computations.Defaults to the maximum number of threads available.
#' @param save Logical. If TRUE, the dataframe of correlations for the dataset is saved as a .csv file.
#' @param filename filename for saved csv file. Only used if save=TRUE. If not specified then the dataset object name is used.
#' @return A dataframe of correlation values. The column genes represent the original clusters whilst the rows represent lagged clusters.
#' @examples
#'
#'
#' @export

coranalysisclusterdataset <- function(cluster.dataset, lag = 0, nthreads = NULL, save = TRUE, filename = NULL) {
    if (is.null(filename) == TRUE) {
        filename <- deparse(substitute(cluster.dataset))  # If a filename isn't specified then the name of the dataframe object is used
    }
    
    
    library(foreach)
    clusters <- unique(cluster.dataset$cluster)  # Vector of cluster numbers
    loading.values <- CircadianTools::loading_gen(length(clusters))  # Calculate milestones for printing progress
    
    correlationdf <- foreach(i = 1:length(clusters), .combine = cbind) %do% {
        CircadianTools::loading_print(iteration = i, loading.values)  # print progress if surpassed a significant milestone
        temp.df <- CircadianTools::coranalysiscluster(i, cluster.dataset, lag = lag, nthreads = nthreads)
        temp.df[, 2]
    }
    
    rownames(correlationdf) <- clusters  # Give the columns the genenames
    colnames(correlationdf) <- clusters  # Give the rows the genenames
    
    if (save == TRUE) {
        write.csv(correlationdf, paste(filename, ".csv", sep = ""))  # Write as a csv file if save=TRUE
    }
    
    return(correlationdf)
    
    
}
