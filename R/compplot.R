#' compplot:
#' @description Plots two genes from a gene activity dataset
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param gene1 The name of the first gene to be plotted. Must be a string.
#' @param gene2 The name of the second gene to be plotted. Must be a string.
#' @param save Logical. If TRUE, saves plots to working directory. Defaults to FALSE.
#'
#'
#' @return Returns or saves ggplot2 object(s).
#' @examples
#' cor_results <- coranalysis(LauraSingleMap, 'comp100002_c0_seq2')

compplot <- function(gene1, gene2, dataset, save = FALSE) {
    timevector <- maketimevector(dataset)
    newgenematrix <- subset(dataset, sample == gene1)
    newgenematrix <- t(newgenematrix[-1])
    
    
    newgenematrix2 <- subset(dataset, sample == gene2)
    newgenematrix2 <- t(newgenematrix2[-1])
    
    
    combined <- cbind(timevector, newgenematrix, newgenematrix2)
    combined <- data.frame(combined)
    names(combined) <- c("timevector", "gene1_activity", "gene2_activity")
    
    
    mean.list <- rep(0, length(unique(timevector)))
    count <- 1
    for (i in unique(timevector)) {
        bob <- subset(combined, timevector == i, select = gene1_activity)
        mean.list[[count]] <- (mean(bob$gene1_activity))
        count = count + 1
    }
    
    mean.list2 <- rep(0, length(unique(timevector)))
    count <- 1
    for (i in unique(timevector)) {
        bob <- subset(combined, timevector == i, select = gene2_activity)
        mean.list2[[count]] <- (mean(bob$gene2_activity))
        count = count + 1
    }
    
    new_mean_list <- rep(mean.list, as.numeric(table(timevector)))
    new_mean_list <- data.matrix(new_mean_list)
    
    new_mean_list2 <- rep(mean.list2, as.numeric(table(timevector)))
    new_mean_list2 <- data.matrix(new_mean_list2)
    
    combined <- cbind(combined$timevector, combined$gene1_activity, combined$gene2_activity, 
        new_mean_list, new_mean_list2)
    combined <- data.frame(combined)
    names(combined) <- c("timevector", "gene1_activity", "gene2_activity", "mean1", "mean2")
    
    
    graphic <- ggplot2::ggplot(data = combined, ggplot2::aes(timevector)) + ggplot2::geom_point(ggplot2::aes(y = gene1_activity, 
        colour = "first gene"), size = 3, alpha = 0.5) + ggplot2::scale_x_continuous(breaks = unique(timevector)) + 
        ggplot2::geom_line(ggplot2::aes(y = mean1, colour = "first gene"))
    graphic <- graphic + ggplot2::xlab("Time (hours)") + ggplot2::ylab("Trancripts Per Million (TPM)") + 
        ggplot2::geom_point(ggplot2::aes(y = gene2_activity, colour = "second gene"), 
            size = 3, alpha = 0.5) + ggplot2::theme_bw() + ggplot2::geom_line(ggplot2::aes(y = mean2, 
        color = "second gene")) + ggplot2::scale_colour_discrete(name = "Gene", breaks = c("first gene", 
        "second gene"), labels = c(gene1, gene2)) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))
    if (save == TRUE) {
        ggplot2::ggsave(paste(gene1, gene2, "comp.png"), graphic, width = 10, height = 4.5, 
            units = "in")
    }
    
    return(graphic)
}
