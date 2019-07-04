#' RainAnalysis :
#' @description Carries out RAIN analysis on a gene dataset
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @return A dataframe object detailing the result of the rain analysis.
#' @examples
#' results <- RainAnalysis(Laurasmappings)
#'
#' @export

RainAnalysis <- function(dataset, period) {
    dataset <- CircadianTools::GeneClean(dataset)
    timevector <- CircadianTools::MakeTimevector(dataset)
    measure <- as.numeric(table(timevector))
    genenames <- dataset[1]
    dataset <- dataset[-1]
    deltat <- unique(timevector)[2] - unique(timevector)[1]
    results <- rain::rain(t(dataset), deltat = deltat, period = period, measure.sequence = measure)
    results <- cbind(genenames, results)
    return(results)
}
