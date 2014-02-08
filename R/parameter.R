##' @export
"parameter<-" <- function(x,...,value) UseMethod("parameter<-")

##' @S3method parameter<- lvmfit
"parameter<-.lvmfit" <- function(x,...,value) {
  parameter(Model(x),...) <- value
  return(x)
}

##' @S3method parameter<- lvm
"parameter<-.lvm" <- function(x,...,value) {
  if (class(value)[1]=="formula") {
    parameter(x,...) <- all.vars(value)
    return(x)
  }
  ## latent(x,silent=TRUE) <- value
  ## covfix(x,value,NULL) <- 1
  ## intfix(x, value) <- value
  newpar <- rep(0,length(value)); names(newpar) <- value
  newfix <- as.list(value); names(newfix) <- value
  x$expar <- c(x$expar,newpar)
  x$exfix <- c(x$exfix,newfix)
  index(x) <- reindex(x)
  x$attributes$parameter[value] <- TRUE
  return(x)
}

##' @export
parameter <- function(x,...)
    names(unlist(x$attributes$parameter))
