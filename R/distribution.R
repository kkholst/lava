##' @export
"distribution<-" <- function(x,...,value) UseMethod("distribution<-")

##' @export
"distribution" <- function(x,...,value) UseMethod("distribution")

##' @S3method distribution<- lvm
"distribution<-.lvm" <- function(x,variable,...,value) {
  if (class(variable)[1]=="formula")
    variable <- all.vars(variable)
  if (length(variable)==1) {
    var. <- paste("\"", variable, "\"", sep="")
    if (is.numeric(value)) value <- list(value)
    mytext <- paste("c(", paste(paste(var.,"=",expression(value),sep=""),collapse=","),")")
    nodeRenderInfo(Graph(x))$"distribution" <- eval(parse(text=mytext))
    return(x)
  }    
  if ((length(value)!=length(variable) & length(value)!=1))
    stop("Wrong number of values")
  for (i in 1:length(variable))
    if (length(value)==1) {
      distribution(x,variable[i],...) <- value
    } else {
      distribution(x,variable[i],...) <- value[[i]]
    }
  return(x)
  
}

##' @S3method distribution lvm
"distribution.lvm" <- function(x,var,...) {
  nodeRenderInfo(Graph(x))$distribution[var]
}


