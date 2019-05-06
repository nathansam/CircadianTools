#' Turningplot:
#' Fits a spline to a given gene in a given dataset. Finds the turning points. Plots the turning points and spline.
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param genename The name of a gene intended for plotting. Must be a string.
#' @param timelag Lags the time. Usually desired if wanting to start from t=0. Defaults to 0.
#' @param save Logical. If TRUE, saves plot to working directory. Defaults to FALSE.
#' @param print Logical. If TRUE renders plot in the plot viewer. Defaults to TRUE
#' @return Prints or saves ggplot2 object(s)
#' @examples
#' turningplot('comp100000_c0_seq2',LauraSingleMap)


turningplot <- function(genename, dataset, timelag = 0, print = TRUE, save = FALSE) {
    genematrix <- subset(dataset, dataset[1] == genename)
    timevector <- maketimevector(genematrix) - timelag  # makes timevector
    genematrix <- genematrix[-1]  #remove gene name from subset


    genesplinefunc <- splinefun(timevector, genematrix)
    x <- seq(timevector[1], tail(timevector, n = 1), by = 0.001)
    y <- genesplinefunc(x, deriv = 0)  # deriv=0 gives the spline itself
    z <- genesplinefunc(x, deriv = 1)  # This is the first derivative of the spline

    spline.df <- data.frame(x, y, z)

    turning.points <- rootSolve::uniroot.all(genesplinefunc, interval = c(6,
        27), deriv = 1)  # this outputs a vector of the turning points, i.e. the points where the derivative=0

    splineplot <- ggplot2::ggplot(ggplot2::aes(x = x, y = y), data = spline.df) +
        ggplot2::geom_line(color = "#39A5AE", size = 2) + ggplot2::theme_bw() +
        ggplot2::xlab("Time(hours)") + ggplot2::ylab("Transcripts Per Million (TPM)") +
        ggplot2::ggtitle(paste("Gene=", genename, ", Mean difference=", round(mean(diff(turning.points)),
            2), " SD=", round(sd(diff(turning.points)), 2))) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))

    for (i in 1:length(turning.points)) {
        splineplot <- splineplot + ggplot2::geom_vline(xintercept = turning.points[i],
            size = 1.5, alpha = 0.75)
    }

    splineplot <- splineplot + ggplot2::theme(text = ggplot2::element_text(size = 12))



    if (save == TRUE) {
        ggplot2::ggsave(paste("Turning Point of ", genename, ".png"), splineplot,
            width = 10, height = 4.5, units = "in")
    }

    if (print == TRUE) {
        return(splineplot)
    }


}
