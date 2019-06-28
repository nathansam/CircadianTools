#' cosinoranalysispar:
#' @description Parallel Implementation of \link{cosinoranalysis}. Fits cosinor models to transcriptomics data and plots the best-fitting models using ggplot2.
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param period The period of rhythmicity which is being tested for. Defaults to 24 (circadian).
#' @param timelag The value of time for which observations begin. Defaults to 6.
#' @param threshold Set significance value for which genes are considered significant and thus plotted. Defaults to 6e-07
#' @param save Logical. If TRUE, saves plots to working directory. Defaults to FALSE.
#' @param print Logical. If TRUE renders significant genes in the plot viewer. Defaults to TRUE
#' @param df Logical. If TRUE a dataframe containing the results of the cosinor analysis will be returned. Defaults to TRUE.
#' @return Prints or saves ggplot2 object(s). Optionally returns dataframe containing gene name and p values from F-test ranking of cosinor models
#' @examples
#' cosinor_results <- cosinoranalysis(LauraSingleMap)
#'
#' @export


cosinoranalysispar <- function(dataset, period = 24, nthreads = NULL, timelag = 6) {
    if (is.null(nthreads) == TRUE) {
        nthreads <- parallel::detectCores()
    }
    library(foreach)
    if (is.null(nthreads) == TRUE) {
        nthreads <- parallel::detectCores()
    }


    dataset <- geneclean(dataset)
    genenumber <- nrow(dataset)  #number of genes in the dataset
    pvalues <- rep(0, genenumber)  #init list of pvalues
    cosinor.pvalue.df <- data.frame(sample = dplyr::select(dataset, 1), pVal = pvalues)  #first column gene name, second:pvalue
    timevector <- maketimevector(dataset)
    loading_values <- loading_gen(genenumber)



    cl <- parallel::makeForkCluster(nthreads)
    doParallel::registerDoParallel(cl)
    cosinor.pvalue.df <- foreach::foreach(i = 1:genenumber, .combine = rbind) %dopar%
        {
            loading_print(i, loading_values)

            genematrix <- dplyr::filter(dataset, dplyr::row_number() == i)
            sample <- genematrix[1, 1]
            genematrix <- genematrix[-1]
            genematrix <- t(genematrix)
            geneexpression <- data.frame(timevector - timelag, genematrix)
            names(geneexpression) <- c("timevector", "activity")
            cosinormodel <- cosinor::cosinor.lm(activity ~ time(timevector), period = period,
                data = geneexpression)
            pVal <- cosinor2::cosinor.detect(cosinormodel)[4]
            data.frame(sample, pVal)
        }

    return(cosinor.pvalue.df)

}
