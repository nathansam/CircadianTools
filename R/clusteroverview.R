#' clusteroverview :
#' @description Plots the mean and error bars for all clusters across time
#'
#' @param cluster.dataset A transcriptomics dataset where the final column details the cluster the gene belongs to. First column should be gene names. All remaining columns should be expression levels.
#' @param nthreads Number of processor threads for the process. If not specifed then the maximum number of logical cores are used.
#' @param save Logical. If TRUE, saves plots. Defaults to FALSE.
#' @param print Logical. If TRUE renders significant genes in the plot viewer. Defaults to TRUE
#' @param path The directory to be used for saving plots to. Creates the directory if it doesn't already exist. Defaults to cluster_overview
#' @return Prints or saves ggplot2 object(s).
#' @examples
#'
#' @export

clusteroverview <- function(cluster.dataset, nthreads = NULL, print = TRUE, save = FALSE, path = "cluster_overview") {
    if (save == TRUE) {
        if (dir.exists(path) == FALSE) {
            # If save==TRUE then create directory for saved plots if it doesn't already exist
            dir.create(path)
        }
    }
    
    for (i in unique(cluster.dataset$cluster)) {
        plot <- CircadianTools::clusterplot(clusterno = i, cluster.dataset, nthreads = nthreads, print = print, 
            save = save, path = path)
        if (print == TRUE) {
            print(plot)
        }
    }
}
