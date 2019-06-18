#' coranalysis:
#' Ranks correlation between a given gene and all over genes in a dataset. Plots both the given gene and highly correlated genes for a given correlation value
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param genename the name of a gene intended for comparison with all other genes in the dataset. Must be a string.
#' @param threshold Set correlation threshold value for which genes are considered significant and thus plotted. Defaults to 0.9
#' @param save Logical. If TRUE, saves plots to working directory. Defaults to FALSE.
#' @param print Logical. If TRUE renders highly correlated genes in the plot viewer. Defaults to TRUE
#' @param df Logical. If TRUE a dataframe containing the correlations of the given gene with all genes in the dataset is returned. Defaults to TRUE.
#' @return Prints or saves ggplot2 object(s). Optionally returns dataframe containing gene names and correlation values
#' @examples
#' cor_results <- coranalysis(LauraSingleMap, 'comp100002_c0_seq2')


coranalysis <- function(genename, dataset, threshold = 0.9, save = FALSE, print = TRUE, df = TRUE) {
    dataset <- geneclean(dataset) # Remove any rows which shows no gene activity
    genenumber <- nrow(dataset)  # Number of genes in the dataset
    cor.df <- data.frame(sample = dplyr::select(dataset, 1), corvalues=rep(0, genenumber))  #first column gene name, second column correlation value
    timevector <- maketimevector(dataset) # Create vector of time values
    loading_values <- loading_gen(genenumber) # Used for the loading bar
    selectedgene <- activity_select(genename, dataset)
    selectedgenedf <- data.frame(timevector, selectedgene)
    names(selectedgenedf) <- c("timevector", "activity")


    selectedmean.list <- rep(0, length((unique(timevector))))
    count <- 1
    for (i in timevector) {
        genesub <- subset(selectedgenedf, timevector == i, select = activity)
        selectedmean.list[[count]] <- (mean(genesub$activity))
        count = count + 1
    }


    for (i in 1:genenumber) {
        loading_print(i, loading_values)

        genematrix <- dplyr::filter(dataset, dplyr::row_number() == i)
        compgenename <- paste(dataset[i,1])
        genematrix<- activity_select(i, dataset)


        selectedgenedf <- data.frame(timevector, genematrix)

        names(selectedgenedf) <- c("timevector", "activity")


        compmean.list <- rep(0, length((unique(timevector))))
        count <- 1
        for (j in timevector) {
            compgenesub <- subset(selectedgenedf, timevector == j, select = activity)
            compmean.list[[count]] <- (mean(compgenesub$activity))
            count = count + 1
        }

        correlation <- cor(selectedmean.list, compmean.list)
        cor.df[i, 2] <- correlation

        if (correlation > threshold) {
            if (correlation != 1) {

                selectedmeans <- data.frame(timevector, means = selectedmean.list, Gene = rep(genename, length(timevector)))
                compmeans <- data.frame(timevector, means = compmean.list, Gene = rep(compgenename, length(timevector)))
                means.df <- rbind(selectedmeans, compmeans)

                selectedactivity <- data.frame(timevector, activity = selectedgene, Gene = rep(genename, length(timevector)))
                names(selectedactivity)[2] <- paste("activity")
                compactivity <- data.frame(timevector, activity = genematrix, Gene = rep(compgenename, length(timevector)))
                activity.df <- rbind(selectedactivity, compactivity)


                corplot <- ggplot2::ggplot(ggplot2::aes(x = timevector, y = means, color = Gene), data = means.df) + ggplot2::geom_line(size = 1) +
                  ggplot2::geom_point(ggplot2::aes(x = timevector, y = activity), data = activity.df, alpha = 0.5, size = 3) + ggplot2::theme_bw() +
                  ggplot2::ggtitle(paste("Cor=", correlation)) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1)) + ggplot2::scale_color_manual(values = c("black",
                  "#39A5AE")) + ggplot2::theme(text = ggplot2::element_text(size = 12)) + ggplot2::xlab("Time (hours)") + ggplot2::ylab("Trancripts Per Million (TPM)")


                if (save == TRUE) {
                  ggplot2::ggsave(paste("Cor_", genename, "_", compgenename, ".png"), corplot, width = 10, height = 4.5, units = "in")
                }
                if (print == TRUE) {
                  print(corplot)
                }
            }
        }
    }
    if (df == TRUE) {
        return(cor.df)
    }
}
