#' clusterspread:
#' @description Shows how many genes are in each cluster after clustering has been applied.
#' @param cluster.dataset A transcriptomics dataset where the final column details the cluster the gene belongs to. First column should be gene names. All remaining columns should be expression levels.
#' @return A dataframe object. The first column is the cluster number. Second column is how many genes belong to that cluster.
#' @examples
#' clusterstats<-clusterspread(clusterdf)
#'
#' @export

clusterspread <- function(cluster.dataset) {
    library(foreach)
    clustersize.df <- foreach(i = unique(cluster.dataset$cluster), .combine = rbind) %do% {
        # Cluster by cluster
        cluster.subset <- subset(cluster.dataset, cluster == i)  # Find all genes in the ith cluster
        cluster.size <- nrow(cluster.subset)  # Find how many genes in the cluster
        data.frame(cluster = i, cluster.size)  # Create row of dataframe object
    }
    return(clustersize.df)
}
