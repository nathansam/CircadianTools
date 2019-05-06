#' cosinoranalysis:
#' Fits cosinor models to transcriptomics data and plots the best-fitting models using ggplot2.
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param period The period of rhythmicity which is being tested for. Defaults to 24 (circadian).
#' @param timelag The value of time for which observations begin. Defaults to 6.
#' @param threshold Set significance value for which genes are considered significant and thus plotted. Defaults to 6e-07
#' @param save Logical. If TRUE, saves plots to working directory. Defaults to FALSE.
#' @param print Logical. If TRUE renders significant genes in the plot viewer. Defaults to TRUE
#' @param pvalues Logical. If TRUE a dataframe containing the results of the cosinor analysis will be returned. Defaults to TRUE.
#' @return Prints or saves ggplot2 object(s). Optionally returns dataframe containing gene name and p values from F-test ranking of cosinor models
#' @examples
#' cosinoranalysis(LauraSingleMap)


cosinoranalysis <- function(dataset, period = 24, timelag = 6, threshold = 6e-07, 
    save = FALSE, print = TRUE, pvalues = TRUE) {
    dataset <- geneclean(dataset)
    genenumber <- nrow(dataset)  #number of genes in the dataset
    pvalues <- rep(0, genenumber)  #init list of pvalues
    cosinor.pvalue.df <- data.frame(sample = dplyr::select(dataset, 1), pvalues)  #first column gene name, second:pvalue
    timevector <- maketimevector(dataset)
    
    twentyfivepercent <- round(0.25 * genenumber)
    fiftypercent <- round(0.5 * genenumber)
    seventyfivepercent <- round(0.75 * genenumber)
    
    
    for (i in 1:genenumber) {
        if (i == twentyfivepercent) {
            print("Progress at 25%")
        }
        if (i == fiftypercent) {
            print("Progress at 50%")
        }
        if (i == seventyfivepercent) {
            print("Progress at 75%")
        }
        
        genematrix <- dplyr::filter(dataset, dplyr::row_number() == i)
        
        genename <- genematrix[1, 1]
        genematrix <- genematrix[-1]
        genematrix <- t(genematrix)
        geneexpression <- data.frame(timevector - timelag, genematrix)
        names(geneexpression) <- c("timevector", "activity")
        cosinormodel <- cosinor::cosinor.lm(activity ~ time(timevector), period = period, 
            data = geneexpression)
        cosinor.pvalue.df[i, 2] <- cosinor2::cosinor.detect(cosinormodel)[4]
        if (cosinor2::cosinor.detect(cosinormodel)[4] < threshold) {
            cosinorplot <- ggplot.cosinor.lm(cosinormodel, endtime = tail(timevector, 
                n = 1) - timelag) + ggplot2::geom_point(ggplot2::aes(y = activity, 
                x = timevector), data = geneexpression, size = 3, alpha = 0.5, 
                color = "#39A5AE") + ggplot2::ggtitle(paste("Gene=", genename, 
                ", P-value=", round(cosinor2::cosinor.detect(cosinormodel)[4], 
                  10))) + ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1)) + 
                ggplot2::theme(text = ggplot2::element_text(size = 12)) + ggplot2::xlab("Time (hours)") + 
                ggplot2::ylab("Trancripts Per Million (TPM)")
            
            
            if (save == TRUE) {
                ggplot2::ggsave(paste("Gene=", genename, ".png"), cosinorplot, 
                  width = 10, height = 4.5, units = "in")
            }
            if (print == TRUE) {
                print(cosinorplot)
            }
        }
    }
    if (pvalues == TRUE) {
        return(cosinor.pvalue.df)
    }
}
