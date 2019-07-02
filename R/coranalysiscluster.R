#' coranalysiscluster
#' @description Correlates the average activity of a cluster with the average activity of every other cluster.
#' @param cluster.no A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param cluster.dataset A transcriptomics dataset where the final column details the cluster the gene belongs to. First column should be gene names. All remaining columns should be expression levels.
#' @param lag Setting any value other than 0 allows a cluster to be correlated with lagged clusters in the dataset. The number denotes the number of timesteps to lag by.
#' @param nthreads The number of threads to be used for parallel computations. Defaults to the maximum number of threads available.
#' @examples
#'
#'
#' @export

coranalysiscluster <- function(cluster.no, cluster.dataset, lag = 0, nthreads = NULL) {
    library(foreach)
    
    main.time.profile <- CircadianTools::clustertimeprofile(cluster.no, cluster.dataset, nthreads = nthreads)  # Get time profile for the cluster specified
    
    if (lag > 0) {
        main.time.profile <- tail(main.time.profile, n = length(main.time.profile) - lag)  # Lag if required
    }
    if (lag < 0) {
        main.time.profile <- head(main.time.profile, n = length(main.time.profile) - lag)  # Lag if required
    }
    
    cluster.quantity <- max(cluster.dataset$cluster)  # Number of clusters
    correlation.df <- data.frame(seq(1, cluster.quantity), rep(0, cluster.quantity))
    colnames(correlation.df) <- c("cluster", "correlation")
    
    for (j in 1:cluster.quantity) {
        comp.time.profile <- CircadianTools::clustertimeprofile(j, cluster.dataset, nthreads = nthreads)  # Get time profile of cluster being compared with the main cluster
        
        if (lag > 0) {
            comp.time.profile <- head(comp.time.profile, n = length(comp.time.profile) - lag)  # Lag if required
        }
        if (lag < 0) {
            comp.time.profile <- tail(comp.time.profile, n = length(comp.time.profile) - lag)  # Lag if required
        }
        
        compcor <- cor(main.time.profile, comp.time.profile)  # Calculate correlation
        correlation.df[j, 2] <- compcor  # Add correlation to dataframe
    }
    return(correlation.df)
}
