#' CytoscapeFile:
#' @description Converts a correlation dataframe object into a format suitable for cytoscape and saves as a csv file.
#' @param cor.dataset A NxN datframe of correlation values created by \link{coranalysisdataset} or \link{coranalysisclusterdataset}
#' @param nthreads Number of processor threads for the process. If not specifed then the maximum number of logical cores are used.
#' @return Dataframe object in the new format.
#' @export

CytoscapeFile <- function(cor.dataset,filename=NULL,nthreads=NULL){

  if (is.null(filename)==TRUE){
    filename <- deparse(substitute(cor.dataset))
    filename <- paste(filename, "_cytoscape", sep="")
  }

  filename <- paste(filename,".csv",sep="")
  CircadianTools::FileConflict(filename)
  # Checks if a file which will be created already exists. Asks the user if this file should be overwritten.


  library(foreach)  #Required for parallelism

  if (is.null(nthreads) == TRUE) {
    # Set the threads to maximum if none is specified
    nthreads <- parallel::detectCores()
  }

  require(foreach)
  sourcelist <- colnames(cor.dataset) # The columns will act as source nodes in Cytoscape
  targetlist <- rownames(cor.dataset)  # The rows will act as target nodes in Cytoscape


  cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
  doParallel::registerDoParallel(cl)

  finaldf <- foreach(j=1:length(sourcelist), .combine = rbind) %dopar% {

    subsetdf <- cor.dataset[,j] # Select column (correlation values for the source)

    sourcedf <- foreach(i=1:length(targetlist), .combine=rbind) %do% {
      data.frame(sourcelist[j], targetlist[i], subsetdf[i]) # Save source, target gene names and correlation values as a row.
      # All rows are combined in returned dataframe
    }
  }
  colnames(finaldf)<-c("source", "target", "correlation")
  parallel::stopCluster(cl)
  write.csv(finaldf, filename, row.names = FALSE) #Write to csv
  return(finaldf)
}
