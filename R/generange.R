#' generange
#' @description Finds the range of gene activity for each gene in a dataframe. The median for the replicates is used for each time point.
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param nthreads Number of processor threads for the filtering. If not specifed then the maximum number of logical cores are used.
#' @examples
#' rangedf <- generange(Laurasmappings, nthreads=4)

generange <- function(dataset, nthreads = NULL) {

    library(foreach)  #Required for parallelism

    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }

    timevector <- CircadianTools::maketimevector(dataset) # List of time values (repeated for replicates)

    genenumber <- nrow(dataset)

    cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
    doParallel::registerDoParallel(cl)


    rangedf <- foreach(i = 1:genenumber, .combine = rbind) %dopar% {
        # Parallel for loop to create a dataframe of gene names and their respective ranges
        gene <- dplyr::filter(dataset, dplyr::row_number() == i)  # Get gene by row
        genename <- gene[, 1]
        genematrix <- t(gene[-1])  # Gene activity as column
        activity.df <- data.frame(genematrix, timevector)
        colnames(activity.df) <- c("activity", "timevector")
        medlist <- c()  #List of medians for the gene
        for (i in unique(timevector)) {
            genesubset <- subset(activity.df, timevector == i) # Populate the median list
            medlist <- c(medlist, median(genesubset$activity))
        }
        generange <- max(medlist) - min(medlist) # Get range
        data.frame(genename, generange)
    }
    parallel::stopCluster(cl)
    return(rangedf)

}
