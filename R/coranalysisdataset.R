#' coranalysisdataset:
#'
#' @description Correlates every gene in a dataset with every other gene in the same dataset. Allows a timelag between genes to be correlated.
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param average The average to be used for comparing the time points. Either "median" or "mean".
#' @param lag Setting any value other than 0 allows a gene to be correlated with lagged genes in the dataset. The number denotes the number of timesteps to lag by.
#' @param par Logical. If TRUE, the correlation for each gene will be generated using multiple threads(in parallel).
#' @param nthreads The number of threads to be used for parallel computations. Only used if par=TRUE. Defaults to the maximum number of threads available.
#' @param save Logical. If TRUE, the dataframe of correlations for the dataset is saved as a .csv file.
#' @param filename filename for saved csv file. Only used if save=TRUE. If not specified then the dataset object name is used.
#' @return A dataframe of correlation values. The column genes represent the original genes whilst the rows represent lagged genes.
#' @examples
#' subdf<-tfilter(Laurasmappings)
#' coranalysisdataset <- coranalysis(subdf, lag=1,filename="cor_tfiltered")
#'
#' @export
coranalysisdataset<-function(dataset, average="median",lag=0, par=TRUE, nthreads=NULL, save=TRUE, filename=NULL){

if (is.null(filename)==TRUE){
  filename<-deparse(substitute(dataset)) # If a filename isn't specified then the name of the dataframe object is used
}

  library(foreach)
  genenames<-as.vector(dataset[,1]) # Vector of names for every gene
  loading.values<-loading_gen(length(genenames)) # Calculate milestones for printing progress

  correlationdf<-foreach(i=1:length(genenames), .combine=cbind) %do% {
    loading_print(iteration=i, loading.values) # print progress if surpassed a significant milestone

    # Calculate correlation for the ith gene with all ( possibly lagged) genes.
    # First column is the name of the gene compared with whilst the second column is the correlations
    if (par==TRUE){ results<-coranalysispar(genename=genenames[i],
              dataset=dataset,lag=lag, average = average, nthreads = nthreads)} # Run in parallel if TRUE
    if (par==FALSE){results<-coranalysis(genename=genenames[i],
                                         dataset=dataset,lag=lag, average = average)} # Don't run in parallel if FALSE
    data.frame(results[,2]) # Add the list of correlations as a column to the correlationdf
  }

  rownames(correlationdf)<-genenames # Give the columns the genenames
  colnames(correlationdf)<-genenames # Give the rows the genenames
  if (save==TRUE){
    write.csv(correlationdf, paste(filename,".csv", sep = "")) # Write as a csv file if save=TRUE
  }
  return(correlationdf) # Return the correlation dataframe.
}
