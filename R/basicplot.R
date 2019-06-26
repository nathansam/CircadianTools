#' basicplot:
#' @description Plots activity data as points and average activity as lines
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param genename The name of a gene intended for plotting. Must be a string.
#' @param timelag Shifts the plot to earlier in time.
#' @param method How should the average activity for each time point be calculated? "median" or "mean"
#' @param save Logical. If TRUE, saves plot to working directory. Defaults to FALSE.
#' @param print Logical. If TRUE renders plot in the plot viewer. Defaults to TRUE
#' @param path The directory to be used for saving plots to. Uses the working directory by default. Not used if save=FALSE

#' @return Prints or saves ggplot2 object(s)
#' @examples
#' basicplot('comp101252_c0_seq2',LauraSingleMap)

basicplot <- function(genename, dataset, timelag = 0, method="median", print = TRUE, save=FALSE, path=NULL) {
    if (save == TRUE){
        if(is.null(path)==FALSE){
        if (dir.exists(path) == FALSE) {
            # If save==TRUE then create directory for saved plots if it doesn't already
            # exist
            dir.create(path)
        }
    }
    }

    genematrix <- subset(dataset, dataset[1] == genename)  # Select gene
    timevector <- (maketimevector(genematrix) - timelag)  # Makes vector of time values
    genematrix <- t(genematrix[-1])  #remove gene name

    genedata <- data.frame(timevector, genematrix)
    names(genedata) <- c("timevector", "activity")


    average.list <- rep(0, length((unique(timevector))))  #Make a list to hold the averaged values for gene activity at each time point.
    count <- 1
    for (i in timevector) {
        genesub <- subset(genedata, timevector == i, select = activity)
        if (method=="mean"){average.list[[count]] <- (mean(genesub$activity))} #compute the mean of the time points
        if (method=="median"){average.list[[count]] <- (median(genesub$activity))} #compute the median of the time points
        count = count + 1
    }

    genedata <- cbind(genedata$timevector, genedata$activity, average.list)
    genedata <- data.frame(genedata)
    names(genedata) <- c("timevector", "activity", "new.average")
    graphic <- ggplot2::ggplot(data = genedata, ggplot2::aes(x = timevector, y = activity)) +
        ggplot2::geom_line(ggplot2::aes(x = timevector, y = new.average), size = 1, color="#412d6b") + ggplot2::geom_point(size = 3,
        alpha = 0.5, color = "#008dd5")
    graphic <- graphic + ggplot2::xlab("Time (Hours)") + ggplot2::ylab("Transcripts Per Million (TPM)") +
        ggplot2::theme_bw()
    graphic <- graphic + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1)) +
        ggplot2::theme(text = ggplot2::element_text(size = 12))
    graphic <- graphic + ggplot2::ggtitle(paste("Gene = ", genename))

    if (save == TRUE) {
            ggplot2::ggsave(paste(genename, ".png"), graphic,path=path ,width = 10, height = 4.5, units = "in")

        }

    if (print == TRUE) {
        return(graphic)
    }
}

