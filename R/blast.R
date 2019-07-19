#' ContigGen
#' @description Finds all unique contig IDs in a transcriptomics dataset
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#'   All other columns should be expression levels.
#' @return A character vector of contig IDs
#' @examples
#' contigs <- ContigGen(Laurasmappings)
#'
#' @export
ContigGen <- function (dataset){
  contigs <- dataset$sample
  contigs <- gsub("_.*","",contigs)
  contigs <-unique(contigs)
  return(contigs)
}
