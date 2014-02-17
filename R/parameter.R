##' @export
"parameter<-" <- function(x,...,value) UseMethod("parameter<-")

##' @S3method parameter<- lvmfit
"parameter<-.lvmfit" <- function(x,...,value) {
  parameter(Model(x),...) <- value
  return(x)
}


##' @S3method parameter<- lvm
"parameter<-.lvm" <- function(x,constrain,start,...,value) {
  if (class(value)[1]=="formula") value <- all.vars(value)
  if (!missing(start)) {
      if (length(start) != length(value)) stop("'start' and 'value' should be of the same lengths")
      start <- as.list(start)
      names(start) <- value
  } else {
      start <- as.list(rep(0,length(value))); names(start) <- value
  }
  if (!missing(constrain)) {      
      newfix <- constrain
      if (!is.list(newfix)) newfix <- as.list(newfix)
  } else {
      newfix <- as.list(value);
  }
  names(newfix) <- value
  x$expar[value] <- start
  x$exfix[value] <- newfix
  index(x) <- reindex(x)
  x$attributes$parameter[value] <- TRUE
  return(x)
}

##' @export
parameter <- function(x,var,...) {
    if (missing(var)) return (names(unlist(x$attributes$parameter)))
    parameter(x,...) <- var
}
