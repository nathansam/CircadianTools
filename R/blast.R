#' ContigGen
#' @description Finds all unique contig IDs in a transcriptomics dataset
#' @param dataset A transcriptomics dataset. First columns should be gene names.
#' All other columns should be expression levels.
#' @return A character vector of contig IDs
#' @examples
#' contigs <- ContigGen(Laurasmappings)
#'
#' @export
ContigGen <- function (dataset){
  contigs <- dataset$sample # Get sample IDs
  # Split strings by underscores into columns
  contigs <- read.table(text = contigs, sep = "_")

  # Join the characters before the first underscore and after the last
  contigs <- paste (contigs[,1], contigs[,3], sep = "_")

  contigs <-unique(contigs)
  return(contigs)
}
