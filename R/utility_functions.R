#' FileConflict
#' @description Checks if a file which will be created already exists and, if necessary asks the user if this file should be overwritten.
#' @param filename The name of a file (with extension included)
#' @export

FileConflict <- function(filename){

  if (file.exists(filename) == TRUE) { # If the file already exists, should it be overwritten?
    prompt <- cat(paste("The file, ", filename, ", already exists. Do you wish to overwrite the file? [Y/n]  \n"))
    continue <- readline(prompt = prompt) # Ask the user

    if (continue == "y" | continue == "Y" | continue == "yes" | continue == "Yes") {
      cat("Overwriting File")
      file.remove(filename) # Overwrite the file if yes
    } else {
      cat("Command Terminated")
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt)) # Stop the function
      stop()
    }
  }
}



#' MedList
#' @description Provides a dataframe of median values at each time point for each gene from a transcriptomics dataset.
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param nthreads Number of processor threads for the filtering. If not specifed then the maximum number of logical cores are used.
#' @examples
#' mediandf <- MedList(Laurasmappings, nthreads=4)
#'
#' @export

MedList <- function(dataset, nthreads = NULL) {

  library(foreach)  #Required for parallelism
  timevector <- CircadianTools::MakeTimevector(dataset)  # List of time values (repeated for replicates)

  genenumber <- nrow(dataset)

  if (is.null(nthreads) == TRUE) {
    # Set the threads to maximum if none is specified
    nthreads <- parallel::detectCores()
  }

  cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
  doParallel::registerDoParallel(cl)

  medlistdf <- foreach(i = 1:genenumber, .combine = rbind) %dopar% {
    # Parallel for loop to create a dataframe of gene names and their respective ranges
    gene <- dplyr::filter(dataset, dplyr::row_number() == i)  # Get gene by row
    genename <- gene[, 1]
    genematrix <- t(gene[-1])  # Gene activity as column
    activity.df <- data.frame(genematrix, timevector)
    colnames(activity.df) <- c("activity", "timevector")
    med.list <- c()  #List of medians for the gene
    for (j in unique(timevector)) {
      genesubset <- subset(activity.df, timevector == j)  # Populate the median list
      med.list <- c(med.list, median(genesubset$activity))
    }
    t(data.frame(med.list))

  }
  parallel::stopCluster(cl)
  colnames(medlistdf) <- unique(timevector)  # Columns of the returned dataframe is the time point
  rownames(medlistdf) <- dataset[, 1]  #Row name is gene name
  return(medlistdf)
}




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





#' Plot a cosinor model
#'
#' Adapted from the Cosinor package by Michael Sachs. Given a cosinor.lm model fit, generate a plot of the data with the fitted values.
#' Optionally allows for plotting by covariate levels 0 and 1. Unlike the original version, the function will plot for the full timecourse rather just one full period.
#'
#'
#' @param object An object of class \code{cosinor.lm}
#' @param x_str Character vector naming the covariate(s) to be plotted. May be NULL to plot overall curve
#' @param endtime The last time value for the time course
#'
#'
#' @export
#'
#'
#'
ggplot.cosinor.lm <- function(object, x_str = NULL, endtime) {

  timeax <- seq(0, endtime, length.out = 200)
  covars <- grep("(rrr|sss)", attr(object$fit$terms, "term.labels"), invert = TRUE, value = TRUE)

  newdata <- data.frame(time = timeax, rrr = cos(2 * pi * timeax/object$period), sss = sin(2 * pi * timeax/object$period))
  for (j in covars) {
    newdata[, j] <- 0
  }
  if (!is.null(x_str)) {

    for (d in x_str) {

      tdat <- newdata
      tdat[, d] <- 1
      newdata <- rbind(newdata, tdat)

    }
    newdata$levels <- ""
    for (d in x_str) {

      newdata$levels <- paste(newdata$levels, paste(d, "=", newdata[, d]))

    }


  }


  newdata$Y.hat <- predict(object$fit, newdata = newdata)

  if (missing(x_str) || is.null(x_str)) {

    ggplot2::ggplot(newdata, ggplot2::aes_string(x = "time", y = "Y.hat")) + ggplot2::geom_line(size = 1.2,
                                                                                                color = "#ffa630")

  } else {

    ggplot2:ggplot(newdata, aes_string(x = "time", y = "Y.hat", col = "levels")) + ggplot2::geom_line()

  }
}


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



#' TAnalysis:
#' @description  A t-test is carried out on gene activity levels between time points and the number of significant increases & decreases is returned.
#' @param row.no The row number of the gene of interest
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param psignificance The maximum p-value for which a result of a t-test is classed as significant
#' @return The number of signficant increases and decreases as a vector
#' @examples
#' ups.downs <- TAnalysis(row.no = 1, dataset = Laurasmappings,  psignificance = 0.01)
#'
#' @export

TAnalysis <- function(row.no, dataset, psignificance = 0.05) {
  timevector <- CircadianTools::MakeTimevector(dataset)  # vector of timevalues
  unique.timevector <- unique(timevector)  #vector of unique timevalues (ignores repetitions)
  gene <- dplyr::filter(dataset, dplyr::row_number() == row.no)  # select gene by row number
  gene <- gene[-1]  #remove sample column
  gene <- as.matrix(t(gene))  # transpose gene to column
  genedf <- data.frame(gene, timevector)  # create dataframe for gene activity and timevector
  names(genedf) <- c("activity", "timevector")
  n.ups <- 0  # Set variable to list the number of significant positive changes found via the t-tests.
  n.downs <- 0  # Set variable to list the number of significant negative changes found via the t-tests.

  for (i in 1:(length(unique.timevector) - 1)) {
    group1 <- subset(genedf, timevector == unique.timevector[i])
    group2 <- subset(genedf, timevector == unique.timevector[i + 1])

    zerogroup <- FALSE  # Flag for if both groups entirely consists of zero values. Assume false
    if (all(group1[, 1] == 0)) {
      # Is group 1 entirely zero values?  Is group 2 entirely zero values?
      if (all(group2[, 1] == 0)) {
        zerogroup <- TRUE  # Both groups are entirely zero values. Don't attempt t-test
      }
    }

    if (zerogroup == FALSE) {
      t.test.obj <- t.test(group1[, 1], group2[, 1], var.equal = TRUE)  # Carry out t-test assuming equal variance
      if (t.test.obj$p.value < psignificance) {
        # If p-value from t-test is below significance threshold Is the t statistic positive?
        if (t.test.obj$statistic > 0) {
          n.ups <- n.ups + 1
        } else {
          # Must be negative then
          n.downs <- n.downs + 1
        }
      }
    }
  }

  return(c(n.ups, n.downs))  # Return the number of signficant increases and decreases
}


#' GeneClean:
#' @description Removes columns and rows which show no gene activity over time
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @return A 'cleaned' transcriptomics dataset
#' @examples
#' cleaned.df<-GeneClean(Laurasmappings)
#'
#' @export


GeneClean <- function(dataset) {
  dataset <- dataset[rowSums(dataset[, -1]) > 0, ]  #removes row(s) which has 0 gene activity
  # dataset <- dataset[colSums(dataset[, -1]) > 0, ] #removes column(s) which has 0 gene activity
  rownames(dataset) <- seq(1, length(dataset[, 1]))
  return(dataset)
}

#' GeneRange
#' @description Finds the range of gene activity for each gene in a dataframe. The median for the replicates is used for each time point.
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param nthreads Number of processor threads for the filtering. If not specifed then the maximum number of logical cores are used.
#' @examples
#' rangedf <- GeneRange(Laurasmappings, nthreads=4)
#'
#' @export

GeneRange <- function(dataset, nthreads = NULL) {

  if (is.null(nthreads) == TRUE) {
    # Set the threads to maximum if none is specified
    nthreads <- parallel::detectCores()
  }

  mediandf <- CircadianTools::MedList(dataset = dataset, nthreads = nthreads)
  genenumber <- nrow(mediandf)

  cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
  doParallel::registerDoParallel(cl)

  rangedf <- foreach(i = 1:genenumber, .combine = rbind) %dopar% {
    gene <- dplyr::filter(mediandf, dplyr::row_number() == i)  # Get gene by row
    generange <- max(gene) - min(gene)
    genename <- dataset[i, 1]
    data.frame(genename, generange)
  }

  parallel::stopCluster(cl)
  return(rangedf)
}


#' GeneScale:
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


#' GeneSub
#' @description Takes an object where the first column is genenames (IE a column of known Circadian genes) and subsets from a dataset containing activity for these genes
#' @param subdf An object where the first column is gene names of interes (IE a column of known Circadian genes or genes found to be signicant)
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param nthreads Number of processor threads used. If not specifed then the maximum number of logical cores are used.
#' @examples
#' newdf <- GeneSub(circadian, Laurasmappings, nthreads=4)
#'
#' @export

GeneSub <- function(subdf, dataframe, nthreads = NULL) {

  library(foreach)  #Required for parallelism

  if (is.null(nthreads) == TRUE) {
    # Set the threads to maximum if none is specified
    nthreads <- parallel::detectCores()
  }

  cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
  doParallel::registerDoParallel(cl)


  newdf <- foreach(i = 1:nrow(subdf), .combine = rbind) %dopar% {
    subset(dataframe, sample == paste(subdf[i, 1]))  # For each gene name in subdf, pull from activity dataframe
  }
  parallel::stopCluster(cl)
  return(newdf)
}


