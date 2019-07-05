#' SingletonNameFinder
#' @description Finds the genes which belong to singleton clusters.
#' @param cluster.dataset A transcriptomics dataset where the final column details the cluster the gene belongs to. First column should be gene names. All remaining columns should be expression levels.
#' @return Returns a vector of gene names for all genes belonging to singleton clusters in the provided dataset.
#' @examples
#' filter.df <- CombiFilter(Laurasmappings)
#' pam.df <- PamClustering(filter.df, 40)
#' singletons <- SingletonNameFinder(pam.df)
#'
#'
#' @export
SingletonNameFinder <- function(cluster.dataset){
  singleton.genes <- c() # Initialise list of genes in singleton clusters
  for (i in unique(cluster.dataset$cluster)){
    sub.df <- subset(cluster.dataset, cluster==i) # Get cluster
    if (nrow(sub.df)==1){ # If only one gene in a cluster
      singleton.genes <- c(singleton.genes, sub.df[1,1]) # Add genename
    }
  }
  return(singleton.genes)
}


#' CommonSingletonFinder
#' @description Finds the genes which belong to common singleton clusters in two different clustered datasets.
#' @param cluster.dataset1 A transcriptomics dataset where the final column details the cluster the gene belongs to. First column should be gene names. All remaining columns should be expression levels.
#' @param cluster.dataset2 A transcriptomics dataset in the same format as cluster.dataset1 but generated via a different clustering method or with different parameters.
#' @examples
#' filter.df <- CombiFilter(Laurasmappings)
#' pam.df <- PamClustering(filter.df, 40)
#' hclust.df <- HClustering(filter.df, 40)
#' common.singletons <- CommonSingletonFinder(pam.df, hclust.df)
#'
#' @export
CommonSingletonFinder <- function(cluster.dataset1, cluster.dataset2){
  singleton.genes1 <- SingletonNameFinder(cluster.dataset1) # Get list of genes in singleton clusters
  singleton.genes2 <- SingletonNameFinder(cluster.dataset2)
  common.genes <- intersect(singleton.genes1, singleton.genes2) # Find the genes which are found in both lists

  cat(paste("There are ", length(common.genes), " singleton genes found in both datasets. \n"))
  cat(paste("There are ", length(singleton.genes1) - length(common.genes), " unique singleton genes in each dataset. \n"))

  return(common.genes)

}
