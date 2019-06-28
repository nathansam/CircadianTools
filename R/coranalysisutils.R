#' @export
activity_select <- function(ID, dataset) {

    if (typeof(ID) == "character") {
        gene <- subset(dataset, dataset[1] == ID)
    } else {
        gene <- dplyr::filter(dataset, dplyr::row_number() == ID)
    }
    gene <- gene[-1]
    gene <- t(gene)
    gene <- matrix(gene)
    return(gene)
}

