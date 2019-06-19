#' geneclean:
#' @description Removes columns and rows which show no gene activity over time
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @return A 'cleaned' transcriptomics dataset
#' @examples
#' geneclean(LauraSingleMap)


geneclean <- function(dataset) {
    dataset <- dataset[rowSums(dataset[, -1]) > 0, ]  #removes row(s) which has 0 gene activity
    dataset <- dataset[colSums(dataset[, -1]) > 0, ]  #removes column(s) which has 0 gene activity
    rownames(dataset) <- seq(1, length(dataset[, 1]))
    return(dataset)
}
