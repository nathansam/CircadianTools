#' genesub
#' @description Takes an object where the first column is genenames (IE a column of known Circadian genes) and subsets from a dataset containing activity for these genes
#' @param subdf An object where the first column is gene names of interes (IE a column of known Circadian genes or genes found to be signicant)
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param nthreads Number of processor threads used. If not specifed then the maximum number of logical cores are used.
#' @examples
#' newdf <- genesub(circadian, Laurasmappings, nthreads=4)
#'
#' @export

genesub <- function(subdf, dataframe, nthreads = NULL) {

    library(foreach)  #Required for parallelism

    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }

    cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
    doParallel::registerDoParallel(cl)


    newdf <- foreach(i = 1:nrow(subdf), .combine = rbind) %dopar% {
        subset(dataframe, sample == paste(subdf[i, 1]))  # For each gene name in subdf, pull from activity dataframe
    }
    parallel::stopCluster(cl)
    return(newdf)
}
