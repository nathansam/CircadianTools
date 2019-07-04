#'#' GeneScale:
#' @description Centers/scales every gene in a transcriptomics dataset
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @return A transcriptomics dataset with centered/scaled genes
#' @examples
#' GeneScale(LaurasMappings)
#' @export

GeneScale <- function(dataset, scale = FALSE, center = TRUE) {
    dataset[-1] <- t(scale(t(dataset[-1]), center = center, scale = scale))
    return(dataset)
}
