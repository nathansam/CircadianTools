#' rainanalysis :
#'
#'
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @return A dataframe object detailing the result of the rain analysis.
#' @examples
#' results <- rainanalysis(Laurasmappings)

rainanalysis <- function(dataset) {
    measure <- c(4, 4, 4, 4, 3, 4, 4, 4)
    dataset <- geneclean(dataset)
    genenames <- dataset[1]
    dataset <- dataset[-1]
    results <- rain::rain(t(dataset), deltat = 3, period = 24, measure.sequence = measure)
    results <- cbind(genenames, results)
    return(results)
}
