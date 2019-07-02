#' combifilter
#' @description Filters a transcriptomics dataset by using zerofilter, anovafilter and sizefilter
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param non_zero_num The minimum number of non-zero activity readings for genes in the filtered dataset
#' @param threshold Set the p-value threshold for the filtering
#' @param cutoff Proportion of total number of genes (which show at least some activity) to be removed
#' @param nthreads Number of processor threads for the filtering. If not specifed then the maximum number of logical cores are used.
#' @param anovafilter Logical. If anovafiltering should be used. Defaults to TRUE.
#' @param zerofilter Logical. If zerofiltering should be used. Defaults to TRUE.
#' @param sizefilter Logical. If sizefiltering should be used. Defaults to TRUE.
#' @examples
#' Laurasmappings_filtered <- combifilter(Laurasmappings, nthreads=4)
#'
#' @export

combifilter <- function(dataset, non_zero_num = 4, threshold = 0.05, cutoff = 0.1, anovafilter = TRUE, zerofilter = TRUE, 
    sizefilter = TRUE, nthreads = NULL) {
    
    if (anovafilter == TRUE) {
        dataset <- CircadianTools::anovafilter(dataset = dataset, threshold = threshold, nthreads = nthreads)  # Filter via ANOVA
        invisible(gc())  # Garbage collection to return RAM to OS before next function
    }
    
    if (zerofilter == TRUE) {
        dataset <- CircadianTools::zerofilter(dataset = dataset, non_zero_num = non_zero_num, nthreads = nthreads)  # Remove genes showing very little activity
        invisible(gc())  # Garbage collection to return RAM to OS before next function
    }
    
    if (sizefilter == TRUE) {
        dataset <- CircadianTools::sizefilter(dataset = dataset, cutoff = cutoff, nthreads = nthreads)  # Filter by removing genes with the smallest range
        invisible(gc())  # Garbage collection to return RAM to OS before next function
    }
    return(dataset)
}
