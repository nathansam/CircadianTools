#' @export
DistanceGen <- function (dataset,metric ,nthreads){
  dataset<-CircadianTools::MedList(dataset, nthreads=nthreads)# Calculate the medians at each timepoint

  if (metric =="abs.correlation"){
  distance <- AbsCorDist(dataset)
} else{
  medians.dataset <- data.frame(medians.dataset)
  distance <- parallelDist::parDist(as.matrix(medians.dataset[-1]), method = metric, threads = nthreads)  #Calculate the distance matrix
}
return(distance)
}
