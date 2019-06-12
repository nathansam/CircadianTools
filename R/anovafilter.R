#' anovafilter:
#' filters a genedata via
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#'
#' @examples
#' Laurasmappings_filtered <- cosinoranalysis(Laurasmappings)


anovafilter <- function(dataset){
  genenumber<-nrow(dataset)
  dataset<-geneclean(dataset)
  count<-0

  timevector<-maketimevector(dataset)
  library(foreach)
  doMC::registerDoMC(4)
  test3<-foreach::foreach(i=1:genenumber, .combine = rbind) %dopar%{
      genematrix <- dplyr::filter(dataset, dplyr::row_number() == i)
      genematrix<-t(genematrix[-1])

      test<-data.frame(genematrix, timevector)
      names(test)<-c("genematrix","timevector")

      #print(genematrix)
      #print(timevector)

    tempaov<-aov(lm(genematrix~timevector, data=test))
    pvalue<-summary(tempaov)[[1]][1,5]
    if (pvalue<0.05){
      dplyr::filter(dataset, dplyr::row_number() == i)
}
  }
  return(test3)
}
