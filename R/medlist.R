#' medlist
#' @description Provides a dataframe of median values at each time point for each gene from a transcriptomics dataset.
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param nthreads Number of processor threads for the filtering. If not specifed then the maximum number of logical cores are used.
#' @examples
#' mediandf <- medlist(Laurasmappings, nthreads=4)
#'
#' @export

medlist <- function(dataset, nthreads = NULL) {

    library(foreach)  #Required for parallelism
    timevector <- CircadianTools::maketimevector(dataset)  # List of time values (repeated for replicates)

    genenumber <- nrow(dataset)

    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }

    cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
    doParallel::registerDoParallel(cl)

    medlistdf <- foreach(i = 1:genenumber, .combine = rbind) %dopar% {
        # Parallel for loop to create a dataframe of gene names and their respective ranges
        gene <- dplyr::filter(dataset, dplyr::row_number() == i)  # Get gene by row
        genename <- gene[, 1]
        genematrix <- t(gene[-1])  # Gene activity as column
        activity.df <- data.frame(genematrix, timevector)
        colnames(activity.df) <- c("activity", "timevector")
        med.list <- c()  #List of medians for the gene
        for (j in unique(timevector)) {
            genesubset <- subset(activity.df, timevector == j)  # Populate the median list
            med.list <- c(med.list, median(genesubset$activity))
        }
        med.list

    }
    parallel::stopCluster(cl)
    medlistdf <- data.frame(medlistdf)
    return(medlistdf)
}
