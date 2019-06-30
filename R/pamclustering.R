#' pamclustering:
#' @description Applies PAM (Partitioning around Medoids) clustering to a transcriptomics dataset and appends a cluster column to this dataset for all genes.
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param k The total number of clusters.
#' @param metric The distance metric to be used to calculate the distances between genes. See parallelDist::parDist for all accepted arguments.
#' @param nthreads Number of processor threads to be used for calculating the distance matrix. If not specifed then the maximum number of logical cores are used.
#' @param scale If the gene activity should be scaled before clustering.
#' @param center If the gene activity should be centered before clustering.
#' @return Returns transcriptomics dataset provided with additional cluster column appended denoted which cluster each gene belongs to.
#' @examples
#' pam.df <- pamclustering(Laurasmappings,20)
#' @export

pamclustering <- function(dataset, k, metric = "euclidean", nthreads = NULL, scale = FALSE, center = TRUE) {

    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }
    medians.dataset <- CircadianTools::medlist(dataset, nthreads=nthreads) # Calculate the medians at each timepoint
    medians.dataset[-1] <- scale(medians.dataset[-1], scale = scale, center = center)  # Scale and center the activity data if TRUE
    medians.dataset<-data.frame(medians.dataset)
    distance <- parallelDist::parDist(as.matrix(medians.dataset[-1]), method = metric, threads = nthreads)  #Calculate the distance matrix
    fit <- cluster::pam(distance, k = k)  # Run the clustering proces
    dataset$cluster <- fit$cluster  # Append the cluster column to the dataset
    return(dataset)
}
