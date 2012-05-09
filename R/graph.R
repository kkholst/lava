##' Extract graph
##' 
##' Extract or replace graph object
##' 
##' 
##' @aliases Graph Graph<-
##' @usage
##' 
##' Graph(x, ...)
##' 
##' Graph(x, ...) <- value
##' 
##' @param x Model object
##' @param value New \code{graphNEL} object
##' @param \dots Additional arguments to be passed to the low level functions
##' @author Klaus K. Holst
##' @seealso \code{\link{Model}}
##' @keywords graphs models
##' @export
##' @examples
##' 
##' m <- lvm(y~x)
##' Graph(m)
##'
##' @export
`Graph` <-
function(x,...) UseMethod("Graph")

##' @S3method Graph lvmfit
##' @S3method Graph lvm
`Graph.lvmfit` <- `Graph.lvm` <-
function(x,add=FALSE,...) {
  if ((is.null(x$graph) || length(x$graph)==0) & add) {
    m <- Model(x)
    return(plot(x,noplot=TRUE))
  }
  else return(x$graph)
}

##' @export
"Graph<-" <- function(x,...,value) UseMethod("Graph<-")

##' @S3method Graph<- lvmfit
"Graph<-.lvmfit" <- function(x,...,value) { x$graph <- value; return(x) }
##' @S3method Graph<- lvm
"Graph<-.lvm" <- function(x,...,value) { x$graph <- value; return(x) }

