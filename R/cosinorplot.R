#' CosinorPlot:
#' @description Fits a cosinor model to a given gene in a given dataset and plots the model
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param genename The name of a gene intended for plotting. Must be a string.
#' @param timelag The value of time for which observations begin. Defaults to 6.
#' @param save Logical. If TRUE, saves plot to working directory. Defaults to FALSE.
#' @param print Logical. If TRUE renders plot in the plot viewer. Defaults to TRUE
#' @return Prints or saves ggplot2 object(s)
#' @examples
#' CosinorPlot('comp100000_c0_seq2', Laurasmappings )
#'
#' @export


CosinorPlot <- function(genename, dataset, timelag = 6, period = 24, print = TRUE, save = FALSE) {
    genematrix <- subset(dataset, dataset[1] == genename)
    timevector <- CircadianTools::MakeTimevector(genematrix)  # makes timevector

    genematrix <- genematrix[-1]  #remove gene name from subset
    genematrix <- t(genematrix)
    geneexpression <- data.frame(timevector - timelag, genematrix)
    names(geneexpression) <- c("timevector", "activity")
    cosinormodel <- cosinor::cosinor.lm(activity ~ time(timevector), period = period, data = geneexpression)
    cosinor.plot.object <- CircadianTools::ggplot.cosinor.lm(cosinormodel, endtime = 21) + ggplot2::geom_point(ggplot2::aes(y = activity,
        x = timevector), data = geneexpression, size = 3, alpha = 0.5, color = "#008dd5") + ggplot2::ggtitle(paste("Gene=",
        genename, ", P-value=", round(cosinor2::cosinor.detect(cosinormodel)[4], 10))) + ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1)) +
        ggplot2::theme(text = ggplot2::element_text(size = 12)) + ggplot2::xlab("Time (hours)") + ggplot2::ylab("Trancripts Per Million (TPM)")

    if (save == TRUE) {
        ggplot2::ggsave(paste("Cosinor_", genename, ".png"), cosinor.plot.object, width = 10, height = 4.5, units = "in")
    }
    if (print == TRUE) {
        return(cosinor.plot.object)
    }
}
