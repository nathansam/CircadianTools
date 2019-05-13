#' basicplot:
#' Plots activity data as points and mean data as lines
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param genename The name of a gene intended for plotting. Must be a string.
#' @param timelag Shifts the plot to earlier in time.
#' @param save Logical. If TRUE, saves plot to working directory. Defaults to FALSE.
#' @param print Logical. If TRUE renders plot in the plot viewer. Defaults to TRUE
#' @return Prints or saves ggplot2 object(s)
#' @examples
#' basicplot('comp101252_c0_seq2',LauraSingleMap)

basicplot <- function(genename, dataset, timelag = 0, save = FALSE, print = TRUE) {
    genematrix <- subset(dataset, dataset[1] == genename)
    timevector <- (maketimevector(genematrix) - timelag)  # makes timevector
    genematrix <- t(genematrix[-1])  #remove gene name from subset
    
    genedata <- data.frame(timevector, genematrix)
    names(genedata) <- c("timevector", "activity")
    
    mean.list <- rep(0, length((unique(timevector))))
    count <- 1
    for (i in timevector) {
        genesub <- subset(genedata, timevector == i, select = activity)
        mean.list[[count]] <- (mean(genesub$activity))
        count = count + 1
    }
    
    
    
    genedata <- cbind(genedata$timevector, genedata$activity, mean.list)
    genedata <- data.frame(genedata)
    names(genedata) <- c("timevector", "activity", "newmean")
    graphic <- ggplot2::ggplot(data = genedata, ggplot2::aes(x = timevector, y = activity)) + ggplot2::geom_line(ggplot2::aes(x = timevector, 
        y = newmean), size = 1) + ggplot2::geom_point(size = 3, alpha = 0.5, color = "#39A5AE")
    graphic <- graphic + ggplot2::xlab("Time (Hours)") + ggplot2::ylab("Transcripts Per Million (TPM)") + ggplot2::theme_bw()
    graphic <- graphic + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1)) + ggplot2::theme(text = ggplot2::element_text(size = 12))
    graphic <- graphic + ggplot2::ggtitle(paste("Gene = ", genename))
    
    if (save == TRUE) {
        ggplot2::ggsave(paste(genename, ".png"), graphic, width = 10, height = 4.5, units = "in")
    }
    if (print == TRUE) {
        return(graphic)
    }
}

