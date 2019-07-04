#' HClustering:
#' @description Applies Hierarchical clustering, clustering to a transcriptomics dataset and appends a cluster column to this dataset for all genes.
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param k The total number of clusters.
#' @param metric The distance metric to be used to calculate the distances between genes. See parallelDist::parDist for all accepted arguments.
#' @param nthreads Number of processor threads to be used for calculating the distance matrix. If not specifed then the maximum number of logical cores are used.
#' @param scale If the gene activity should be scaled before clustering.
#' @param center If the gene activity should be centered before clustering.
#' @return Returns transcriptomics dataset provided with additional cluster column appended denoted which cluster each gene belongs to.
#' @examples
#' pam.df <- Hclustering(Laurasmappings,75)
#'
#' @export
HClustering <- function(dataset, k = 10, metric = "euclidean", nthreads = NULL, scale = FALSE, center = TRUE) {

    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }
    dataset <- CircadianTools::GeneScale(dataset, scale = scale, center = center)  # Center / scale the gene activity for each gene
    medians.dataset <- CircadianTools::MedList(dataset, nthreads = nthreads)  # Calculate the medians at each timepoint
    # medians.dataset[-1] <- scale(medians.dataset[-1], scale = scale, center = center) # Scale/center the data
    medians.dataset <- data.frame(medians.dataset)
    distance <- parallelDist::parDist(as.matrix(medians.dataset[-1]), method = metric, threads = nthreads)  #Calculate the distance matrix
    fit <- hclust(distance)  # Run the clustering process
    clusters <- dendextend::cutree(fit, k = k)  # Cut the dendogram such that there are k clusters
    dataset$cluster <- clusters  # Append the cluster column to the dataset
    return(dataset)
}
