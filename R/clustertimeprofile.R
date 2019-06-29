#' clustertimeprofile
#' @description Provides a dataframe of median values at each time point for each cluster.
#' @param cluster.no A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param cluster.dataset A transcriptomics dataset where the final column details the cluster the gene belongs to. First column should be gene names. All remaining columns should be expression levels.
#' @param nthreads The number of threads to be used for parallel computations. Defaults to the maximum number of threads available.
#' @examples
#'
#'
#' @export

clustertimeprofile<-function(cluster.no, cluster.dataset, nthreads=NULL){
  cluster.sub<-subset(cluster.dataset, cluster==cluster.no) # Subset cluster
  cluster.sub$cluster<-NULL # Remove cluster column
  med<-medlist(cluster.sub, nthreads=nthreads) # Return the median activity label for each timepoint of each gene

  timesteps<-ncol(med) # Number of timesteps

  profile<-rep(0, timesteps)

    for (i in 1:timesteps){
      profile[i]<-mean(med[,i]) # Find the mean activity value for the cluster at each time points
    }
  return(profile)
}
