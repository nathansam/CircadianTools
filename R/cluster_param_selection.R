#' PamParamSelection
#' @description Runs PAM with differing numbers of partitions and returns validation metrics.
#' @param dataset A transcriptomics dataset. Preferably filtered first. First columns should be gene names. All other columns should be expression levels.
#' @param k A numeric vector giving the number of clusters to be evaluated.
#' @param nthreads The number of threads to be used for parallel computations. Warning! Increasing this value will massively increase RAM usage. If NULL then the maximum number of threads available will be used.
#' @examples
#' filter.df <- CombiFilter(Laurasmappings)
#' k.options <- seq(10,100, by=10)
#' pam.validation <- PamParamSelection(filterdf, k=k.options)
#' @export
PamParamSelection <- function(dataset, k=c(2,5,10), nthreads=2){

  if (is.null(nthreads) == TRUE) {
    # Set the threads to maximum if NULL is given as argument for nthreads
    nthreads <- parallel::detectCores()
  }

  `%dopar%` <- foreach::`%dopar%` # Load the dopar binary operator from foreach package
  dataset <- CircadianTools::GeneScale(dataset) # Center each gene
  dataset <- data.frame(CircadianTools::MedList(dataset, nthreads=nthreads)) # Reduce dataset to median at each time point.
  distance <- parallelDist::parDist(as.matrix(dataset[-1])) # Calculate distance using euclidean distance
  cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
  doParallel::registerDoParallel(cl)


  result.df <- foreach::foreach(i = k, .combine = rbind ) %dopar% {
    clust<-cluster::pam(distance, k=i) # Run PAM
    dunn <- clValid::dunn(distance, clust$cluster) # Calculate Dunn index
    connect <- (clValid::connectivity(distance, clust$cluster)) # Calculate connectivity
    silhouette <- mean(cluster::silhouette(clust$clustering, distance)) # Calculate Silhouette width
    data.frame(i,dunn, connect, silhouette, "PAM") # Make row of result.df

  }
  parallel::stopCluster(cl)
  colnames(result.df)<-c("k", "Dunn", "Connectivity", "Silhouette", "Method") # Column headings
  return(result.df)
}


#' HclustParamSelection
#' @description Runs hierarchical clustering with differing numbers of partitions and returns validation metrics.
#' @param dataset A transcriptomics dataset. Preferably filtered first. First columns should be gene names. All other columns should be expression levels.
#' @param k A numeric vector giving the number of clusters to be evaluated.
#' @param nthreads The number of threads to be used for parallel computations.If NULL then the maximum number of threads available will be used.
#' filter.df <- CombiFilter(Laurasmappings)
#' k.options <- seq(10,100, by=10)
#' hclust.validation <- hclustParamSelection(filterdf, k=k.options)
#' @export
 HclustParamSelection <- function(dataset, k=c(2,5,10), nthreads=4){

   if (is.null(nthreads) == TRUE) {
     # Set the threads to maximum if none is specified
     nthreads <- parallel::detectCores()
   }

   `%dopar%` <- foreach::`%dopar%` # Load the dopar binary operator from foreach package
   dataset <- CircadianTools::GeneScale(dataset) # Center each gene
   dataset <- data.frame(CircadianTools::MedList(dataset, nthreads=nthreads)) # Reduce dataset to median at each time point.
   distance <- parallelDist::parDist(as.matrix(dataset[-1])) # Calculate distance using euclidean distance
   clust <- hclust(distance) # Run hclustering

   cl <- parallel::makeForkCluster(nthreads)  # Create cluster for parallelism
   doParallel::registerDoParallel(cl)

   result.df <- foreach::foreach(i = k, .combine = rbind ) %dopar% {
     cluster <- dendextend::cutree(clust, k=i) # Cut tree
     dunn <- clValid::dunn(distance, cluster) # Calculate Dunn index
     connect <- (clValid::connectivity(distance, cluster)) # Calculate connectivity
     silhouette <- mean(cluster::silhouette(cluster, distance)) # Calculate Silhouette width
     data.frame(i,dunn, connect, silhouette ,"Hierarchical") # Make row of result.df

   }
   parallel::stopCluster(cl)
   colnames(result.df)<-c("k", "Dunn", "Connectivity", "Silhouette", "Method") # Column headings
   return(result.df)
}

#' ClusterParamSelection
#' @description Calculates validation metrics for different clustering methods and varying the number of partitions. The validation metrics are plotted.
#' @param dataset A transcriptomics dataset. Preferably filtered first. First columns should be gene names. All other columns should be expression levels.
#' @param k A numeric vector giving the number of clusters to be evaluated.
#' @param method The clustering method(s) to be used. Multiple methods can be considered by passing a vector. Currently accepts 'pam' and 'hclust'
#' @param nthreads The number of threads to be used for parallel computations.If NULL then the maximum number of threads available will be used.
#' @param save.plot Logical. If TRUE then saves the plots generated
#' @param save.df Logical. If TRUE then saves the validation metric results as a csv file.
#' @param path The directory to be used for saving plots and the validation metric results to. Uses the name of the dataset object appended with '_validation' if this argument is not specified
#' @examples
#' filter.df <- CombiFilter(Laurasmappings)
#' k.options <- seq(10,100, by=10)
#' hclust.validation <- ClusterParamSelection(filterdf, k=k.options)
#' @export
ClusterParamSelection <- function(dataset,k=c(2,5,10), method=c('pam', 'hclust'), nthreads=2,save.plot=TRUE,save.df=TRUE, path=NULL ){
  if (is.null(path) == TRUE) {
    path <- deparse(substitute(dataset))  # If a filename isn't specified then the name of the dataframe object is used
    path <- paste(path, '_validation', sep="") # Add _validation to directory

  }


  if (dir.exists(path)==FALSE){
    dir.create(path) # Create directory if it doesn't already exist
  }

  validation.df <- data.frame() # Initialize dataframe to hold validation results

  if ("pam" %in% method ==TRUE){
    pam.results <- CircadianTools::PamParamSelection(dataset=dataset, k=k, nthreads=nthreads)
    validation.df <- rbind(validation.df, pam.results) # Add pam validation results if pam is specified in methods
  }

  if ("hclust" %in% method == TRUE){
    hclust.results <- CircadianTools::HclustParamSelection(dataset = dataset, k=k, nthreads = nthreads)
    validation.df <- rbind(validation.df, hclust.results) # Add hclust validation results if hclust is specified in methods
  }

  colours.vector <- c("#008dd5", "#ba1200", "ffa630", "#840032", "#412d6b") # Vector of colours used in package
  colours.vector <- colours.vector[1:length(method)] # Select number of colours required


  #### Dunn ####
  p <- ggplot2::ggplot(ggplot2::aes(x=k, y=Dunn, color= Method), data=validation.df)+
    ggplot2::geom_line(size=1.2)+ggplot2::theme_bw()+ ggplot2::ggtitle("Dun Index Plot")+
    ggplot2::scale_colour_manual(values=colours.vector) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))
  print (p)

  if (save.plot==TRUE){
    ggplot2::ggsave("dunn_index.png", plot=p, path=path, width = 10
                    , height = 4.5, units = "in") # Save the plot
  }


  #### Connectivity ####
  p <- ggplot2::ggplot(ggplot2::aes(x=k, y=Connectivity, color= Method), data=validation.df)+
    ggplot2::geom_line(size=1.2)+ggplot2::theme_bw()+ ggplot2::ggtitle("Connectivity Plot")+
    ggplot2::scale_colour_manual(values=colours.vector) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))
  print (p)

  if (save.plot==TRUE){
    ggplot2::ggsave("connectivity_plot.png", plot=p, path=path, width = 10
                    , height = 4.5, units = "in") # Save the plot
  }


  #### Silhouette ####
  p <- ggplot2::ggplot(ggplot2::aes(x=k, y=Silhouette, color= Method), data=validation.df)+
    ggplot2::geom_line(size=1.2)+ggplot2::theme_bw()+ ggplot2::ggtitle("Silhouette Plot")+
    ggplot2::scale_colour_manual(values=colours.vector)+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))
  print (p)

  if (save.plot==TRUE){
    ggplot2::ggsave("silhouette_plot.png", plot=p, path=path, width = 10
                    , height = 4.5, units = "in") # Save the plot
  }



  if (save.df==TRUE){
    df.filename<-paste(path, "/results.csv", sep="") # Save the validation results as .csv file
    write.csv(validation.df, file=df.filename, row.names=FALSE)
  }

  return(validation.df)
}
