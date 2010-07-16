"parameter<-" <- function(x,...,value) UseMethod("parameter<-")
"parameter<-.lvm" <- function(x,...,value) {
  if (class(value)[1]=="formula") {
    parameter(x,...) <- all.vars(value)
    return(x)
  }
  latent(x) <- value
  covfix(x,value,NULL) <- 1
  intfix(x, value) <- value
  nodeData(Graph(x), value, attr="parameter") <- TRUE
  return(x)
}

parameter <- function(x,...)
  vars(x)[unlist(lapply(nodeData(Graph(x)), function(z) z$parameter))]
