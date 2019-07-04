#' TFilter:
#' @description Experimental! Applies a filter where a t.test is carried out on gene activity levels between time points. The number of significant changes between time points is found. If there is a sufficient number of significant changes and close to as many positive changes as negative changes then the gene is included in the filtered dataset
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param maxdifference The maximum difference between the number of signifcant positive t statistic and the number of signifcant negative t statistic.
#' @param minchanges The minimum number of significant changes found via t-tests for the gene to be included in the filtered dataset.
#' @param nthreads Number of processor threads to be used for calculating the distance matrix. If not specifed then the maximum number of logical cores are used.
#' @param psignificance The maximum p-value for which a result of a t-test is classed as significant
#' @return Returns a filtered transcriptomics dataset
#' @examples
#' filterdf <- TFilter(Laurasmappings)
#'
#' @export

TFilter <- function(dataset, maxdifference = 1, minchanges = 2, psignificance = 0.05, nthreads = NULL) {
    library(foreach)

    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }

    cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
    doParallel::registerDoParallel(cl)

    dataset[-1] <- scale(dataset[-1], scale = FALSE, center = TRUE)

    filterdf <- foreach(i = 1:nrow(dataset), .combine = rbind) %dopar% {
        ups.downs <- TAnalysis(row.no = i, dataset = dataset, psignificance = psignificance)
        # tanalysis returns two values as a vector. First values represents a positive significant change between time
        # points whilst second value represents a negative significant change.
        observed.difference <- abs(ups.downs[1] - ups.downs[2])  # Finds the difference between the number of positive and negative changes
        total.changes <- ups.downs[1] + ups.downs[2]  # The total number of significant changes

        if (total.changes >= minchanges) {
            # if total changes is above the user set minimum if the difference between signficant 'ups' and 'downs' is below
            # the user set maximum
            if (observed.difference <= maxdifference) {
                data.frame(dataset[i, ])  # Add the gene as part of the filtered dataset
            }
        }
    }

    parallel::stopCluster(cl)
    return(filterdf)
}
