#' sizefilter
#' @description Filters the genes with the smallest range from a transcriptomics dataset
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param cutoff Proportion of total number of genes (which show at least some activity) to be removed
#' @param nthreads Number of processor threads for the filtering. If not specifed then the maximum number of logical cores are used.
#' @examples
#' rangedf <- sizefilter(Laurasmappings, nthreads=4)
#'
#' @export

sizefilter <- function(dataset, cutoff = 0.1, nthreads = NULL) {

    dataset <- CircadianTools::geneclean(dataset)  # Remove genes with no activity
    genenumber <- nrow(dataset)
    timevector <- CircadianTools::maketimevector(dataset)

    rangedf <- generange(dataset = dataset, nthreads = nthreads)  # Calculate range of activity for each gene
    rangedf <- rangedf[order(rangedf$generange, decreasing = TRUE), ]  # Order by biggest range
    rank <- round(genenumber * (1 - cutoff))  # Number of genes which will be returned after filtering
    rangedf <- rangedf[1:rank, ]  # Selects the genes with the highest range
    filterdf <- genesub(rangedf, dataset, nthreads)  # Returns genes found in rangedf with their respective activity readings
    row.names(filterdf) <- seq(1, nrow(filterdf))

    return(filterdf)
}
