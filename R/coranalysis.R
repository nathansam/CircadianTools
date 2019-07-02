#' coranalysis:
#' @description Ranks correlation between a given gene and all other genes in a dataset. Plots both the given gene and highly correlated genes for a given correlation value
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param genename the name of a gene intended for comparison with all other genes in the dataset. Must be a string.
#' @param threshold Set correlation threshold value for which genes are considered significant and thus plotted. Defaults to 0.9
#' @param lag Setting any value other than 0 allows a gene to be correlated with lagged genes in the dataset. The number denotes the number of timesteps to lag by.
#' @param average The average to be used for comparing the time points. Either 'median' or 'mean'.
#' @param save Logical. If TRUE, saves plots to working directory. Defaults to FALSE.
#' @param print Logical. If TRUE renders highly correlated genes in the plot viewer. Defaults to TRUE
#' @param df Logical. If TRUE a dataframe containing the correlations of the given gene with all genes in the dataset is returned. Defaults to TRUE.
#' @return Prints or saves ggplot2 object(s). Optionally returns dataframe containing gene names and correlation values
#' @examples
#' cor_results <- coranalysis(LauraSingleMap, 'comp100002_c0_seq2')
#'
#' @export


coranalysis <- function(genename, dataset, threshold = 0.9, average = "median", lag = 0, save = FALSE, print = TRUE, 
    df = TRUE) {
    dataset <- geneclean(dataset)  # Remove any rows which shows no gene activity
    
    if (save == TRUE) {
        directory <- paste("cor_", genename)
        if (dir.exists(directory) == FALSE) {
            dir.create(directory)
        }
    }
    
    genenumber <- nrow(dataset)  # Number of genes in the dataset
    cor.df <- data.frame(sample = dplyr::select(dataset, 1), corvalues = rep(0, genenumber))  #first column gene name, second column correlation value
    timevector <- CircadianTools::maketimevector(dataset)  # Create vector of time values
    loading_values <- CircadianTools::loading_gen(genenumber)  # Used for the loading bar
    selectedgene <- CircadianTools::activity_select(genename, dataset)
    selectedgenedf <- data.frame(timevector, selectedgene)
    names(selectedgenedf) <- c("timevector", "activity")
    
    
    selectedaverage.list <- rep(0, length((unique(timevector))))
    count <- 1
    for (i in unique(timevector)) {
        genesub <- subset(selectedgenedf, timevector == i, select = activity)
        if (average == "mean") {
            selectedaverage.list[[count]] <- (mean(genesub$activity))
        }
        if (average == "median") {
            selectedaverage.list[[count]] <- (median(genesub$activity))
        }
        count = count + 1
    }
    
    if (lag > 0) {
        selectedaverage.list <- tail(selectedaverage.list, n = length(selectedaverage.list) - lag)
    }
    
    if (lag < 0) {
        selectedaverage.list <- head(selectedaverage.list, n = length(selectedaverage.list) - lag)
    }
    
    for (i in 1:genenumber) {
        loading_print(i, loading_values)
        
        genematrix <- dplyr::filter(dataset, dplyr::row_number() == i)
        compgenename <- paste(dataset[i, 1])
        genematrix <- activity_select(i, dataset)
        
        
        
        selectedgenedf <- data.frame(timevector, genematrix)
        
        names(selectedgenedf) <- c("timevector", "activity")
        
        
        compaverage.list <- rep(0, length((unique(timevector))))
        count <- 1
        for (j in unique(timevector)) {
            compgenesub <- subset(selectedgenedf, timevector == j, select = activity)
            if (average == "mean") {
                compaverage.list[[count]] <- (mean(compgenesub$activity))
            }
            if (average == "median") {
                compaverage.list[[count]] <- (median(compgenesub$activity))
            }
            count = count + 1
        }
        
        if (lag > 0) {
            compaverage.list <- head(compaverage.list, n = length(compaverage.list) - lag)
        }
        if (lag < 0) {
            compaverage.list <- tail(compaverage.list, n = length(compaverage.list) - lag)
        }
        
        
        correlation <- cor(selectedaverage.list, compaverage.list)
        cor.df[i, 2] <- correlation
        if (save == TRUE || print == TRUE) {
            if (correlation > threshold) {
                if (correlation != 1) {
                  
                  
                  myplot <- compplot(as.character(genename), compgenename, dataset)
                  myplot <- myplot + ggplot2::ggtitle(paste("Correlation = ", correlation))
                  
                  if (save == TRUE) {
                    ggplot2::ggsave(paste("Cor_", genename, "_", compgenename, ".png"), myplot, path = directory, 
                      width = 10, height = 4.5, units = "in")
                  }
                  if (print == TRUE) {
                    print(myplot)
                  }
                }
            }
        }
    }
    if (df == TRUE) {
        return(cor.df)
    }
}
