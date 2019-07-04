#' AnovaFilter
#' @description Filters a gene activity dataframe via ANOVA
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param nthreads Number of processor threads for the filtering. If not specifed then the maximum number of logical cores are used.
#' @param threshold Set the p-value threshold for the filtering
#' @examples
#' Laurasmappings_filtered <- anovafilter(Laurasmappings, nthreads=4)
#'
#' @export

AnovaFilter <- function(dataset, threshold = 0.05, nthreads = NULL) {
    library(foreach)  #Required for parallelism

    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }

    dataset <- CircadianTools::GeneClean(dataset)  # Remove genes with no activity
    genenumber <- nrow(dataset)  # Number of genes in the dataset
    timevector <- CircadianTools::MakeTimevector(dataset)  # List of time values (repeated for replicates)

    cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
    doParallel::registerDoParallel(cl)  # Register cluster

    filterdf <- foreach(i = 1:genenumber, .combine = rbind) %dopar% {
        # Parallel for loop to create dataframe of significant genes
        gene <- dplyr::filter(dataset, dplyr::row_number() == i)  # Get gene by row
        genematrix <- t(gene[-1])  # Remove gene name
        tempaov <- aov(lm(as.numeric(genematrix) ~ as.factor(timevector)))  # Fit model and create aov object
        pvalue <- summary(tempaov)[[1]][1, 5]  # Get the p-value from the aov object
        if (pvalue < threshold) {
            gene  # Return   the gene (this gene will be combined with other significant genes found in the for loop via rbind to form the dataframe)
        }
    }
    parallel::stopCluster(cl)

    rownames(filterdf) <- seq(1, nrow(filterdf))  #Rebuild the row names

    return(filterdf)
}
