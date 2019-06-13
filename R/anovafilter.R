#' anovafilter:
#' filters a genedata via
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param nthreads Number of processor threads for the filtering. If not specifed then the maximum logical cores are used.
#' @examples
#' Laurasmappings_filtered <- anovafilter(Laurasmappings, nthreads=4)


anovafilter <- function(dataset, nthreads = NULL) {
    library(foreach)
    if (is.null(nthreads) == TRUE) {
        nthreads <- parallel::detectCores()
    }
    
    dataset <- geneclean(dataset)
    genenumber <- nrow(dataset)
    count <- 0
    
    timevector <- maketimevector(dataset)
    cl <- parallel::makeForkCluster(nthreads)
    doParallel::registerDoParallel(cl)
    
    filterdf <- foreach(i = 1:genenumber, .combine = rbind) %dopar% {
        genematrix <- dplyr::filter(dataset, dplyr::row_number() == i)
        genematrix <- t(genematrix[-1])
        
        test <- data.frame(genematrix, timevector)
        names(test) <- c("genematrix", "timevector")
        tempaov <- aov(lm(genematrix ~ timevector, data = test))
        pvalue <- summary(tempaov)[[1]][1, 5]
        if (pvalue < 0.05) {
            dplyr::filter(dataset, dplyr::row_number() == i)
        }
    }
    return(filterdf)
}
