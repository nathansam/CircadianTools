#' generange
#' @description Finds the range of gene activity for each gene in a dataframe. The median for the replicates is used for each time point.
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param nthreads Number of processor threads for the filtering. If not specifed then the maximum number of logical cores are used.
#' @examples
#' rangedf <- generange(Laurasmappings, nthreads=4)
#'
#' @export

generange <- function(dataset, nthreads = NULL) {
    
    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }
    
    mediandf <- CircadianTools::medlist(dataset = dataset, nthreads = nthreads)
    genenumber <- nrow(mediandf)
    
    cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
    doParallel::registerDoParallel(cl)
    
    rangedf <- foreach(i = 1:genenumber, .combine = rbind) %dopar% {
        gene <- dplyr::filter(mediandf, dplyr::row_number() == i)  # Get gene by row
        generange <- max(gene) - min(gene)
        genename <- dataset[i, 1]
        data.frame(genename, generange)
    }
    
    parallel::stopCluster(cl)
    return(rangedf)
}
