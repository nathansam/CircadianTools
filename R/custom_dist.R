
#' @export

AbsCorDist <- function (dataset, nthreads=NULL){
  if( is.null(nthreads)==TRUE){
    nthreads <- parallel::detectCores()
  }
  fn <- function(x,y) {1-abs(cor(dataset[x,],dataset[y,]))}
  distances <- proxy::dist(1:nrow(dataset), method=fn)
  remove(fn)
return(distances)
}
