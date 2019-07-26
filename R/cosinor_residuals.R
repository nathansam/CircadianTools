#' CosinorResidualPlot
#' @description Fits a cosinor model to a gene and plot the residuals
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#' All other columns should be expression levels.
#' @param genename The name of a gene intended for plotting. Must be a string.
#' @param timelag Shifts the plot to earlier in time.
#' @param period The period of rhythmicity which is being tested for.
#' Defaults to 24 (circadian).
#' @param save Logical. If TRUE, saves plot to working directory. Defaults to
#' FALSE.
#' @param print Logical. If TRUE renders plot in the plot viewer. Defaults to
#' TRUE
#' @param path The directory to be used for saving plots to. Uses the working
#' directory by default. Not used if save=FALSE
#' @return Prints or saves ggplot2 object(s)
#' @examples
#' CosinorResidualPlot('comp99801_c1_seq1', Laurasmappings)
#'
#' @export
CosinorResidualPlot <- function(genename, dataset, timelag = 6, period = 24, print = TRUE, save = FALSE,
    path = NULL) {

    if (save == TRUE) {
        if (is.null(path) == FALSE) {
            if (dir.exists(path) == FALSE) {
                # If save==TRUE then create directory for saved plots if needed
                dir.create(path)
            }
        }
    }

    genematrix <- subset(dataset, dataset[1] == genename)
    timevector <- CircadianTools::MakeTimevector(genematrix)  # makes timevector

    genematrix <- genematrix[-1]  #remove gene name from subset
    genematrix <- t(genematrix)
    geneexpression <- data.frame(timevector - timelag, genematrix)
    names(geneexpression) <- c("timevector", "activity")
    cosinormodel <- cosinor::cosinor.lm(activity ~ time(timevector), period = period, data = geneexpression)
    # plot(timevector,cosinormodel$fit$residuals)

    test <- data.frame(geneexpression, resids = cosinormodel$fit$residuals)
    p <- ggplot2::ggplot(ggplot2::aes(x = timevector, y = resids), data = test) + ggplot2::geom_point(size = 3,
        alpha = 0.5, color = "#008dd5")
    p <- p + ggplot2::xlab("Time (Hours)") + ggplot2::ylab("Residuals") + ggplot2::theme_bw()

    p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1)) + ggplot2::theme(text = ggplot2::element_text(size = 12))
    p <- p + ggplot2::ggtitle(paste("Residual Plot for ", genename, " with ", period, " hour period",
        sep = ""))

    if (save == TRUE) {
        ggplot2::ggsave(paste(genename, ".png"), p, path = path, width = 10, height = 4.5, units = "in")
    }


    if (print == TRUE) {
        return(p)
    }
}
