#' coranalysisdataset:
#'
#' @description Correlates every gene in a dataset with every other gene in the same dataset. Allows a timelag between genes to be correlated.
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param average The average to be used for comparing the time points. Either 'median' or 'mean'.
#' @param lag Setting any value other than 0 allows a gene to be correlated with lagged genes in the dataset. The number denotes the number of timesteps to lag by.
#' @param nthreads The number of threads to be used for parallel computations. Defaults to the maximum number of threads available.
#' @param save Logical. If TRUE, the dataframe of correlations for the dataset is saved as a .csv file.
#' @param filename filename for saved csv file. Only used if save=TRUE. If not specified then the dataset object name is used.
#' @return A dataframe of correlation values. The column genes represent the original genes whilst the rows represent lagged genes.
#' @examples
#' subdf<-tfilter(Laurasmappings)
#' cordf <- coranalysisdataset(subdf, lag=1,filename='cor_tfiltered')
#'
#' @export
coranalysisdataset <- function(dataset, average = "median", lag = 0, nthreads = NULL, save = TRUE, filename = NULL) {
    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }

    if (is.null(filename) == TRUE) {
        # If a filename isn't specified then the name of the dataframe object is used
        filename <- deparse(substitute(dataset))
    }

    library(foreach)
    cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
    doParallel::registerDoParallel(cl)

    genenames <- as.vector(dataset[, 1])  # Vector of names for every gene

    correlationdf <- foreach(i = 1:length(genenames), .combine = cbind) %dopar% {

        # Calculate correlation for the ith gene with all ( possibly lagged) genes.
        results <- coranalysis(genename = genenames[i], dataset = dataset, lag = lag, average = average, print=FALSE)

        data.frame(results[, 2])  # Add the list of correlations as a column to the correlation dataframe
    }

    rownames(correlationdf) <- genenames  # Give the columns the genenames
    colnames(correlationdf) <- genenames  # Give the rows the genenames
    if (save == TRUE) {
        write.csv(correlationdf, paste(filename, ".csv", sep = ""))  # Write as a csv file if save=TRUE
    }
    parallel::stopCluster(cl)
    return(correlationdf)  # Return the correlation dataframe.
}
