#' clusterplot :
#' @description Plots the mean and error bars for the genes in a cluster across time
#'
#' @param clusterno The number which identifies the cluster
#' @param cluster.dataset A transcriptomics dataset where the final column details the cluster the gene belongs to. First column should be gene names. All remaining columns should be expression levels.
#' @param nthreads Number of processor threads for the process. If not specifed then the maximum number of logical cores are used.
#' @param save Logical. If TRUE, saves plots. Defaults to FALSE.
#' @param print Logical. If TRUE renders significant genes in the plot viewer. Defaults to TRUE
#' @param path The directory to be used for saving plots to. Uses the working directory by default. Not used if save=FALSE
#' @return Prints or saves ggplot2 object(s).
#' @examples
#'
#' @export

clusterplot <- function(clusterno, cluster.dataset, nthreads = NULL, print = TRUE,
    save = FALSE, path = NULL) {
    library(foreach)
    subdf <- subset(cluster.dataset, cluster == clusterno)  # Subset by cluster
    subdf$cluster <- NULL  #remove the cluster column

    unique.time.vector <- unique(maketimevector(subdf))  # Get the time values
    subdfmedians <- medlist(subdf, nthreads = nthreads)  # Generates the median at each time point for each gene


    if (nrow(subdfmedians) != 1) {
        single.gene.cluster = FALSE  # Set logical flag for there being more than one gene in the cluster (Error bars and standard deviation is required)
        } else {
        single.gene.cluster = TRUE  # Set logical flag for there being 1 gene in the cluster (Error bars and standard deviation is not required)
    }

    graphdf <- foreach(i = 1:ncol(subdfmedians), .combine = rbind) %do% {

        column <- subdfmedians[,i ]  # Select all values per column (per timepoint)
        meanval <- mean(column)  # Calculate the mean value for this timepoint
        time <- unique.time.vector[i]  # Find the actual value of time for this timepoint

        if (single.gene.cluster == TRUE) {
            data.frame(time, meanval)  # Store just the time value and mean if only one gene
        } else {
            # If more than one gene in the cluster
            se <- sd(column)/length(column)  # Calculate the standard error for this time point

            data.frame(time, meanval, se)  # Store time value, mean and standard deviation
        }
    }

    p <- ggplot2::ggplot(graphdf, ggplot2::aes(x = time, y = meanval))  # Create the ggplot2 object

    if (single.gene.cluster == FALSE) {
        p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = meanval - (2*se), ymax = meanval +
            (2*se)), width = 1.5, size = 1, position = ggplot2::position_dodge(0.05),
            color = "#ba1200", alpha = 0.7)  # Add error bars if more than 1 gene in cluster
    }


    p <- p + ggplot2::geom_line(size = 1, color="#412d6b") + ggplot2::geom_point(size = 4, color = "#008dd5") +
        ggplot2::xlab("Time (Hours)") + ggplot2::ylab("Transcripts Per Million (TPM)") +
        ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1)) +
        ggplot2::theme(text = ggplot2::element_text(size = 12)) + ggplot2::ggtitle(paste("Cluster = ",
        clusterno))  # Add line, points and change appearance to match packaged appearance


    if (save == TRUE) {
        ggplot2::ggsave(paste("cluster=", clusterno, ".png"), p, path = path, width = 10,
            height = 4.5, units = "in")
    }

    if (print == TRUE) {
        return(p)
    }

}
