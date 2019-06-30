#' datasetplot:
#' @description Saves plots of all genes in a dataset. WARNING! Don't run on a large dataset! Intended for a filtered dataset
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param timelag Shifts the plot to earlier in time.
#' @param method How should the average activity for each time point be calculated? "median" or "mean"
#' @param nthreads Number of processor threads to be used for calculating the distance matrix. If not specifed then the maximum number of logical cores are used.
#' @param path The directory to be used for saving plots to. Uses the working directory by default and uses the name of the dataset object if this argument is nots specified)  Not used if save=FALSE
#'
#' @export

datasetplot <- function(dataset, timelag = 0,method="median", nthreads = NULL, path = NULL) {

  if (is.null(path)==TRUE){
    path<-deparse(substitute(dataset)) # If a filename isn't specified then the name of the dataframe object is used
  }

  library(foreach)
   if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }

    cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
    doParallel::registerDoParallel(cl)

    genenames <- dataset[, 1]
    filterdf <- foreach(i = 1:length(genenames)) %dopar% {
        basicplot(genenames[i], dataset = dataset, timelag = timelag,method=method ,print = FALSE, save = TRUE, path = path)
    }



}
