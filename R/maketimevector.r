#' MakeTimevector:
#' @description Produces a vector of time values for the gene activity readings.
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @return A vector of time values for the genes
#' @examples
#' MakeTimevector(Laurasmappings)
#'
#' @export

MakeTimevector <- function(dataset) {
    columnnames <- colnames(dataset)
    columnnames <- columnnames[-1]  #remove first column (sample)
    timevector <- gsub("CT", "", columnnames)  #remove CT characters
    timevector <- gsub("\\..*", "", timevector)  #remove characters after full stop
    timevector <- as.numeric(timevector)  #converts from string to numeric
    return(timevector)
}

