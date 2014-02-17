##' @export
"parameter<-" <- function(x,...,value) UseMethod("parameter<-")

##' @S3method parameter<- lvmfit
"parameter<-.lvmfit" <- function(x,...,value) {
  parameter(Model(x),...) <- value
  return(x)
}

##' @S3method parameter<- lvm
"parameter<-.lvm" <- function(x,start,...,value) {
  if (class(value)[1]=="formula") {
    parameter(x,...) <- all.vars(value)
    return(x)
  }
  ## latent(x,silent=TRUE) <- value
  ## covfix(x,value,NULL) <- 1
  ## intfix(x, value) <- value
  if (!missing(start)) {
      if (length(start) != length(value)) stop("'start' and 'value' should be of the same lengths")
      start <- as.list(start)
      names(start) <- value
  } else {
      start <- as.list(rep(0,length(value))); names(start) <- value
  }
  newfix <- as.list(value); names(newfix) <- value
  x$expar[value] <- start
  x$exfix[value] <- newfix
  ##x$expar <- c(x$expar,start)
  ##x$exfix <- c(x$exfix,newfix)
  index(x) <- reindex(x)
  x$attributes$parameter[value] <- TRUE
  return(x)
}

##' @export
parameter <- function(x,...)
    names(unlist(x$attributes$parameter))
