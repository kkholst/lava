##' Extract graph
##' 
##' Extract or replace graph object
##' 
##' 
##' @aliases Graph Graph.lvmfit Graph.lvm Graph<- Graph<-.lvmfit Graph<-.lvm
##' @usage
##' Graph(x, ...)
##' Graph(x, ...) <- value
##' @param x Model object
##' @param value New \code{graphNEL} object
##' @param \dots Additional arguments to be passed to the low level functions
##' @return Returns a graph object (\code{graphNEL})
##' @author Klaus K. Holst
##' @seealso \code{\link{Model}}
##' @keywords graphs models
##' @examples
##' 
##' m <- lvm(y~x)
##' Graph(m)
##'
##' @export
`Graph` <-
function(x,...) UseMethod("Graph")

`Graph.lvmfit` <- `Graph.lvm` <-
function(x,add=FALSE,...) {
  if ((is.null(x$graph) || length(x$graph)==0) & add) {
    m <- Model(x)
    return(plot(x,noplot=TRUE))
  }
  else return(x$graph)
}

"Graph<-" <- function(x,...,value) UseMethod("Graph<-")
"Graph<-.lvmfit" <- function(x,...,value) { x$graph <- value; return(x) }
"Graph<-.lvm" <- function(x,...,value) { x$graph <- value; return(x) }

