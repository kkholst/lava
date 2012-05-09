##' @export
`survival` <- function(x,...) UseMethod("survival")
##' @export
"survival<-" <- function(x,...,value) UseMethod("survival<-")

##' @S3method survival<- lvm
"survival<-.lvm" <- function(x,...,value) {
  survival(x, value, ...)
}

##' @S3method survival lvm
`survival.lvm` <-
function(x,var=NULL, ...) {
  if (is.null(var)) {
    survidx <- unlist(nodeData(Graph(x), attr="survival"))
    if (length(survidx)>0)
      return(names(survidx)[survidx])
    else
      NULL    
  }
  x <- addattr(x,attr="shape",var=var,val="box")
  nodeData(Graph(x), var, attr="survival") <- TRUE
  nodeData(Graph(x), var, attr="normal") <- TRUE
  return(x)
}

