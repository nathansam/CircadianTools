#' coranalysispar:
#' @description Parallel Implementation of \link{coranalysis}. Ranks correlation between a given gene and all over genes in a dataset. Plots both the given gene and highly correlated genes for a given correlation value
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param genename the name of a gene intended for comparison with all other genes in the dataset. Must be a string.
#' @param nthreads Number of processor threads for the process. If not specifed then the maximum number of logical cores are used.
#' @return Returns dataframe containing gene names and correlation values
#' @examples
#' cor_results <- coranalysis('comp100002_c0_seq2',Laurasmappings)


coranalysispar <- function(genename, dataset, nthreads = NULL) {

    library(foreach)  #Required for parallelism
    if (is.null(nthreads) == TRUE) {
        # Set the threads to maximum if none is specified
        nthreads <- parallel::detectCores()
    }

    dataset <- geneclean(dataset)  # Remove any rows which shows no gene activity
    genenumber <- nrow(dataset)  # Number of genes in the dataset
    cor.df <- data.frame(sample = dplyr::select(dataset, 1), corvalues = rep(0, genenumber))  #first column gene name, second column correlation value
    timevector <- maketimevector(dataset)  # Create vector of time values
    selectedgene <- activity_select(genename, dataset)
    selectedgenedf <- data.frame(timevector, selectedgene)
    names(selectedgenedf) <- c("timevector", "activity")


    selectedmean.list <- rep(0, length((unique(timevector))))
    count <- 1
    for (i in unique(timevector)) {
        genesub <- subset(selectedgenedf, timevector == i, select = activity)
        selectedmean.list[[count]] <- (mean(genesub$activity))
        count = count + 1
    }

    cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
    doParallel::registerDoParallel(cl)

    cor.df <- foreach(i = 1:genenumber, .combine = rbind) %dopar% {

        genematrix <- dplyr::filter(dataset, dplyr::row_number() == i)
        compgenename <- paste(dataset[i, 1])
        genematrix <- activity_select(i, dataset)


        selectedgenedf <- data.frame(timevector, genematrix)

        names(selectedgenedf) <- c("timevector", "activity")


        compmean.list <- rep(0, length((unique(timevector))))
        count <- 1
        for (j in unique(timevector)) {
            compgenesub <- subset(selectedgenedf, timevector == j, select = activity)
            compmean.list[[count]] <- (mean(compgenesub$activity))
            count = count + 1
        }

        correlation <- cor(selectedmean.list, compmean.list)
        data.frame(compgenename, correlation)


    }

    parallel::stopCluster(cl)
    return(cor.df)
}
