
#' @export

AbsCorDist <- function (dataset, nthreads=NULL){
  if( is.null(nthreads)==TRUE){
    nthreads <- parallel::detectCores()
  }

  meds<-MedList(dataset, nthreads=nthreads)
  fn <- function(x,y) {1-abs(cor(meds[x,],meds[y,]))}
  distances <- proxy::dist(1:nrow(dataset), method=fn)
  remove(fn)
return(distances)
}
