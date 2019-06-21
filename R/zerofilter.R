#' zerofilter
#' @description Filters a transcriptomics dataset such that there is a minimum number of non-zero activity readings for each gene.
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param non_zero_num The minimum number of non-zero activity readings for genes in the filtered dataset
#' @param nthreads Number of processor threads for the filtering. If not specifed then the maximum number of logical cores are used.
#' @examples
#' Laurasmappings_filtered <- zerofilter(Laurasmappings, nthreads=4)

zerofilter <- function(dataset, non_zero_num = 4, nthreads = NULL) {
    dataset <- CircadianTools::geneclean(dataset)
    library(foreach)
    if (is.null(nthreads) == TRUE) {
        nthreads <- parallel::detectCores()  # Set the threads to maximum if none is specified
    }
    
    genenumber <- nrow(dataset)
    timevector <- maketimevector(dataset)
    cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
    doParallel::registerDoParallel(cl)
    filterdf <- foreach(i = 1:genenumber, .combine = rbind) %dopar% {
        num_zeroes <- sum(as.list((dataset[i, ][-1])) == 0)  # Number of zeroes
        if (num_zeroes < (length(timevector) - non_zero_num)) {
            dataset[i, ]
        }
    }
    parallel::stopCluster(cl)
    rownames(filterdf) <- seq(1, nrow(filterdf))  #Rebuild the row names
    return(filterdf)
}
