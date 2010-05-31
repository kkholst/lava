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

