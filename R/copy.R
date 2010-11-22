"copy<-" <- function(object,...,value) UseMethod("copy<-")
"copy<-.lvm" <- function(object,...,value) {
  if (class(value)[1]=="formula") {
    value <- all.vars(value)
  }
  if (length(value)<2) stop("Provide name of new variable")
  
  regression(object,value[-1],value[1]) <- 1
  covariance(object,value[-1]) <- 0  
  intercept(object,value[-1]) <- 0
  return(object)
}
