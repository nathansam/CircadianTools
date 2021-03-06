#' AbsCorDist
#' @description Calculates a distance matrix based on the distance measure of:
#' \deqn{1 - |cor(x, y)|}
#' @param dataset A dataset created my using \link{MedList} on a transcriptomics
#'  dataset
#' @param nthreads Number of processor threads for the filtering. If not
#'  specifed then the maximum number of logical cores are used.
#' @examples
#' meds <- MedList(Laurasmappings, nthreads = 2)
#' distance <- AbsCorDist(meds, nthreads = 2)
#' @export

AbsCorDist <- function(dataset, nthreads = NULL) {
    if (is.null(nthreads) == TRUE) {
        nthreads <- parallel::detectCores()
    }
    fn <- function(x, y) {
        1 - abs(stats::cor(dataset[x, ], dataset[y, ]))
    }
    distances <- proxy::dist(1:nrow(dataset), method = fn)
    remove(fn)
    return(distances)
}
