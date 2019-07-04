#' clusters.to.text
#' @description Takes a dataframe of clusters and stores the name of all genes in a text file. The row number deontes the cluster number.
#' @param cluster.dataset A transcriptomics dataset where the final column details the cluster the gene belongs to. First column should be gene names. All remaining columns should be expression levels.
#' @param filename The filename of the saved text file. If not given then the name of the correlation dataframe object will be used. The.txt extension is not needed.
#' @export
clusters.to.text <- function(cluster.dataset, filename = NULL) {

    if (is.null(filename) == TRUE) {
        filename <- deparse(substitute(cluster.dataset)) # If filename is not given then use name of cluster.dataset object
    }

    filename <- paste(filename, ".txt", sep = "") # Add .txt extension to the filename

    file.conflict(filename) # Checks if a file which will be created already exists. Asks the user if this file should be overwritten.


    cluster.list <- unique(cluster.dataset$cluster) # vector of cluster 'names'

    for (i in cluster.list) {
        subdf <- subset(cluster.dataset, cluster == i) # Subset by cluster
        line <- subdf[, 1]  # genenames (first column) # Get the gene names for this cluster
        size <- length(line) # How many genenames for the cluster
        line <- matrix(line, ncol = size) # Create matrix to act as row in the txt file
        write(line, file = filename, append = TRUE, ncolumns = size, sep = ",") # Write the row to the txt file
    }

}
