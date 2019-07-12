#' ClusterDatasetPlot :
#' @description Plots the mean and error bars for all clusters across time
#'
#' @param cluster.dataset A transcriptomics dataset where the final column details the cluster the gene belongs to. First column should be gene names. All remaining columns should be expression levels.
#' @param nthreads Number of processor threads for the process. If not specifed then the maximum number of logical cores are used.
#' @param save Logical. If TRUE, saves plots. Defaults to FALSE.
#' @param print Logical. If TRUE renders significant genes in the plot viewer. Defaults to TRUE
#' @param path The directory to be used for saving plots to. Creates the directory if it doesn't already exist. Defaults to cluster_overview
#' @return Prints or saves ggplot2 object(s).
#' @examples
#' filter.df <- CombiFilter(Laurasmappings)
#' pam.df <- PamClustering(filter.df, k = 75)
#' ClusterDatasetPlot(pam.df)
#'
#' @export

ClusterDatasetPlot <- function(cluster.dataset, nthreads = NULL, print = TRUE, save = FALSE, path = "cluster_overview") {
    if (save == TRUE) {
        if (dir.exists(path) == FALSE) {
            # If save==TRUE then create directory for saved plots if it doesn't already exist
            dir.create(path)
        }
    }
    
    for (i in unique(cluster.dataset$cluster)) {
        plot <- CircadianTools::ClusterPlot(clusterno = i, cluster.dataset, nthreads = nthreads, print = print, 
            save = save, path = path)
        if (print == TRUE) {
            print(plot)
        }
    }
}

#' ClusterPlot :
#' @description Plots the mean and error bars for the genes in a cluster across time
#'
#' @param clusterno The number which identifies the cluster
#' @param cluster.dataset A transcriptomics dataset where the final column details the cluster the gene belongs to. First column should be gene names. All remaining columns should be expression levels.
#' @param nthreads Number of processor threads for the process. If not specifed then the maximum number of logical cores are used.
#' @param save Logical. If TRUE, saves plots. Defaults to FALSE.
#' @param print Logical. If TRUE renders significant genes in the plot viewer. Defaults to TRUE
#' @param path The directory to be used for saving plots to. Uses the working directory by default. Not used if save=FALSE
#' @return Prints or saves a ggplot2 object.
#' @examples
#' filter.df <- CombiFilter(Laurasmappings)
#' pam.df <- PamClustering(filter.df, k = 75)
#' ClusterPlot(2, pam.df)
#' @export

ClusterPlot <- function(clusterno, cluster.dataset, nthreads = NULL, print = TRUE, save = FALSE, path = NULL) {
    `%do%` <- foreach::`%do%`  # Load the do binary operator from foreach package
    subdf <- subset(cluster.dataset, cluster == clusterno)  # Subset by cluster
    subdf$cluster <- NULL  #remove the cluster column
    
    unique.time.vector <- unique(CircadianTools::MakeTimevector(subdf))  # Get the time values
    subdfmedians <- CircadianTools::MedList(subdf, nthreads = nthreads)  # Generates the median at each time point for each gene
    
    if (nrow(subdfmedians) != 1) {
        single.gene.cluster = FALSE  # Set logical flag for there being more than one gene in the cluster (Error bars and standard deviation is required)
    } else {
        single.gene.cluster = TRUE  # Set logical flag for there being 1 gene in the cluster (Error bars and standard deviation is not required)
    }
    
    graphdf <- foreach::foreach(i = 1:ncol(subdfmedians), .combine = rbind) %do% {
        
        column <- subdfmedians[, i]  # Select all values per column (per timepoint)
        meanval <- mean(column)  # Calculate the mean value for this timepoint
        time <- unique.time.vector[i]  # Find the actual value of time for this timepoint
        
        if (single.gene.cluster == TRUE) {
            data.frame(time, meanval)  # Store just the time value and mean if only one gene
        } else {
            # If more than one gene in the cluster
            se <- sd(column)/sqrt(length(column))  # Calculate the standard error for this time point
            
            data.frame(time, meanval, se)  # Store time value, mean and standard deviation
        }
    }
    
    p <- ggplot2::ggplot(graphdf, ggplot2::aes(x = time, y = meanval))  # Create the ggplot2 object
    
    if (single.gene.cluster == FALSE) {
        p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = meanval - (2 * se), ymax = meanval + (2 * se)), 
            width = 1.5, size = 1, position = ggplot2::position_dodge(0.05), color = "#ba1200", alpha = 0.7)  # Add error bars if more than 1 gene in cluster
    }
    
    
    p <- p + ggplot2::geom_line(size = 1, color = "#412d6b") + ggplot2::geom_point(size = 4, color = "#008dd5") + 
        ggplot2::xlab("Time (Hours)") + ggplot2::ylab("Transcripts Per Million (TPM)") + ggplot2::theme_bw() + 
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1)) + ggplot2::theme(text = ggplot2::element_text(size = 12)) + 
        ggplot2::ggtitle(paste("Cluster = ", clusterno))  # Add line, points and change appearance to match packaged appearance
    
    
    if (save == TRUE) {
        ggplot2::ggsave(paste("cluster=", clusterno, ".png"), p, path = path, width = 10, height = 4.5, 
            units = "in")
    }
    
    if (print == TRUE) {
        return(p)
    }
    
}


#' ClusterSpread:
#' @description Shows how many genes are in each cluster after clustering has been applied.
#' @param cluster.dataset A transcriptomics dataset where the final column details the cluster the gene belongs to. First column should be gene names. All remaining columns should be expression levels.
#' @return A dataframe object. The first column is the cluster number. Second column is how many genes belong to that cluster.
#' @examples
#' filter.df <- CombiFilter(Laurasmappings)
#' pam.df <- PamClustering(filter.df, k = 75)
#' clusterstats<-ClusterSpread(pam.df)
#' @export

ClusterSpread <- function(cluster.dataset) {
    `%do%` <- foreach::`%do%`  # Load the do binary operator from foreach package
    clustersize.df <- foreach::foreach(i = unique(cluster.dataset$cluster), .combine = rbind) %do% {
        # Cluster by cluster
        cluster.subset <- subset(cluster.dataset, cluster == i)  # Find all genes in the ith cluster
        cluster.size <- nrow(cluster.subset)  # Find how many genes in the cluster
        data.frame(cluster = i, cluster.size)  # Create row of dataframe object
    }
    return(clustersize.df)
}


#' ClusterText
#' @description Takes a dataframe of clusters and stores the name of all genes in a text file. The row number deontes the cluster number.
#' @param cluster.dataset A transcriptomics dataset where the final column details the cluster the gene belongs to. First column should be gene names. All remaining columns should be expression levels.
#' @param filename The filename of the saved text file. If not given then the name of the correlation dataframe object will be used. The.txt extension is not needed.
#' @examples
#' filter.df <- CombiFilter(Laurasmappings)
#' pam.df <- PamClustering(filter.df, k = 75)
#' ClusterText(pam.df)
#' @export
ClusterText <- function(cluster.dataset, filename = NULL) {
    
    if (is.null(filename) == TRUE) {
        filename <- deparse(substitute(cluster.dataset))  # If filename is not given then use name of cluster.dataset object
    }
    
    filename <- paste(filename, ".txt", sep = "")  # Add .txt extension to the filename
    
    file.conflict(filename)  # Checks if a file which will be created already exists. Asks the user if this file should be overwritten.
    
    
    cluster.list <- unique(cluster.dataset$cluster)  # vector of cluster 'names'
    
    for (i in cluster.list) {
        subdf <- subset(cluster.dataset, cluster == i)  # Subset by cluster
        line <- subdf[, 1]  # genenames (first column) # Get the gene names for this cluster
        size <- length(line)  # How many genenames for the cluster
        line <- matrix(line, ncol = size)  # Create matrix to act as row in the txt file
        write(line, file = filename, append = TRUE, ncolumns = size, sep = ",")  # Write the row to the txt file
    }
}



#' ClusterTimeProfile
#' @description Provides a dataframe of median values at each time point for each cluster.
#' @param cluster.no A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param cluster.dataset A transcriptomics dataset where the final column details the cluster the gene belongs to. First column should be gene names. All remaining columns should be expression levels.
#' @param nthreads The number of threads to be used for parallel computations. Defaults to the maximum number of threads available.
#' @examples
#' filter.df <- CombiFilter(Laurasmappings)
#' pam.df <- PamClustering(filter.df, k = 75)
#' time.profile<-ClusterTimeProfile(1,pam.df)
#'
#' @export

ClusterTimeProfile <- function(cluster.no, cluster.dataset, nthreads = NULL) {
    cluster.sub <- subset(cluster.dataset, cluster == cluster.no)  # Subset cluster
    cluster.sub$cluster <- NULL  # Remove cluster column
    med <- CircadianTools::MedList(cluster.sub, nthreads = nthreads)  # Return the median activity label for each timepoint of each gene
    
    timesteps <- ncol(med)  # Number of timesteps
    
    profile <- rep(0, timesteps)
    
    for (i in 1:timesteps) {
        profile[i] <- mean(med[, i])  # Find the mean activity value for the cluster at each time points
    }
    return(profile)
}


#' CommonSingletonFinder
#' @description Finds the genes which belong to common singleton clusters in two different clustered datasets.
#' @param cluster.dataset1 A transcriptomics dataset where the final column details the cluster the gene belongs to. First column should be gene names. All remaining columns should be expression levels.
#' @param cluster.dataset2 A transcriptomics dataset in the same format as cluster.dataset1 but generated via a different clustering method or with different parameters.
#' @examples
#' filter.df <- CombiFilter(Laurasmappings)
#' pam.df <- PamClustering(filter.df, k = 40)
#' hclust.df <- HClustering(filter.df, k = 40)
#' common.singletons <- CommonSingletonFinder(pam.df, hclust.df)
#'
#' @export
CommonSingletonFinder <- function(cluster.dataset1, cluster.dataset2) {
    singleton.genes1 <- SingletonNameFinder(cluster.dataset1)  # Get list of genes in singleton clusters
    singleton.genes2 <- SingletonNameFinder(cluster.dataset2)
    common.genes <- intersect(singleton.genes1, singleton.genes2)  # Find the genes which are found in both lists
    
    cat(paste("There are ", length(common.genes), " singleton genes found in both datasets. \n"))
    cat(paste("There are ", length(singleton.genes1) - length(common.genes), " unique singleton genes in each dataset. \n"))
    
    return(common.genes)
    
}


#' HClustering:
#' @description Applies Hierarchical clustering to a transcriptomics dataset and appends a cluster column to this dataset for all genes.
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param k The total number of clusters.
#' @param metric The distance metric to be used to calculate the distances between genes. See parallelDist::parDist for all accepted arguments.
#' @param nthreads Number of processor threads to be used for calculating the distance matrix. If not specifed then the maximum number of logical cores are used.
#' @param scale If the gene activity should be scaled before clustering.
#' @param center If the gene activity should be centered before clustering.
#' @return Returns transcriptomics dataset provided with additional cluster column appended denoted which cluster each gene belongs to.
#' @examples
#' pam.df <- Hclustering(Laurasmappings, k = 75)
#'
#' @export
HClustering <- function(dataset, k = 10, metric = "euclidean", nthreads = NULL, scale = FALSE, center = TRUE) {
    
    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }
    dataset <- CircadianTools::GeneScale(dataset, scale = scale, center = center)  # Center / scale the gene activity for each gene
    
    if (metric == "abs.correlation") {
        distance <- AbsCorDist(dataset)
    } else {
        medians.dataset <- CircadianTools::MedList(dataset, nthreads = nthreads)  # Calculate the medians at each timepoint
        distance <- parallelDist::parDist(as.matrix(medians.dataset), method = metric, threads = nthreads)  #Calculate the distance matrix
    }
    
    fit <- hclust(distance)  # Run the clustering process
    clusters <- dendextend::cutree(fit, k = k)  # Cut the dendogram such that there are k clusters
    dataset$cluster <- clusters  # Append the cluster column to the dataset
    return(dataset)
}


#' PamClustering:
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
#' pam.df <- PamClustering(Laurasmappings,k = 20)
#' @export

PamClustering <- function(dataset, k, metric = "euclidean", nthreads = NULL, scale = FALSE, center = TRUE) {
    
    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }
    dataset <- CircadianTools::GeneScale(dataset, scale = scale, center = center)  # Center / scale the gene activity for each gene
    
    if (metric == "abs.correlation") {
        distance <- AbsCorDist(dataset)
    } else {
        medians.dataset <- CircadianTools::MedList(dataset, nthreads = nthreads)  # Calculate the medians at each timepoint
        distance <- parallelDist::parDist(as.matrix(medians.dataset), method = metric, threads = nthreads)  #Calculate the distance matrix
    }
    fit <- cluster::pam(distance, k = k)  # Run the clustering proces
    dataset$cluster <- fit$cluster  # Append the cluster column to the dataset
    return(dataset)
}


#' SingletonNameFinder
#' @description Finds the genes which belong to singleton clusters.
#' @param cluster.dataset A transcriptomics dataset where the final column details the cluster the gene belongs to. First column should be gene names. All remaining columns should be expression levels.
#' @return Returns a vector of gene names for all genes belonging to singleton clusters in the provided dataset.
#' @examples
#' filter.df <- CombiFilter(Laurasmappings)
#' pam.df <- PamClustering(filter.df, k = 40)
#' singletons <- SingletonNameFinder(pam.df)
#'
#'
#' @export
SingletonNameFinder <- function(cluster.dataset) {
    singleton.genes <- c()  # Initialise list of genes in singleton clusters
    for (i in unique(cluster.dataset$cluster)) {
        sub.df <- subset(cluster.dataset, cluster == i)  # Get cluster
        if (nrow(sub.df) == 1) {
            # If only one gene in a cluster
            singleton.genes <- c(singleton.genes, sub.df[1, 1])  # Add genename
        }
    }
    return(singleton.genes)
}
