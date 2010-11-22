
`heavytail` <- function(x,...) UseMethod("heavytail")
"heavytail<-" <- function(x,...,value) UseMethod("heavytail<-")

"heavytail<-.lvm" <- function(x,...,value) {
  if (class(value)[1]=="formula") {
    return(heavytail(x,all.vars(value),...))
  }
  binary(x, value, ...)
}

`heavytail.lvm` <-
function(x,var=NULL,df=1,...) {
  if (is.null(var)) {
    htidx <- unlist(nodeData(Graph(x), attr="heavytail"))
    if (length(htidx)>0 && any(htidx!=0)) {
      return(htidx[htidx>0])      
    }
    return(NULL)
  }
  nodeData(Graph(x), var, attr="heavytail") <- df
  return(x)
}


heavytail.init.hook <- function(x,...) {
  nodeDataDefaults(x,"heavytail") <- 0; x
}
heavytail.sim.hook <- function(x,data,...) {
  hvar <- heavytail(x)
  n <- nrow(data)
  count <- 0
  for (i in names(hvar)) {
    df <- hvar[i]
    Z <- rchisq(n, df=df)
    Y <- sqrt(df)*data[,i]/sqrt(Z)
    data[,i] <- Y
##    data[,i] <- data[,i]*sqrt(df/rchisq(n,df=df))
  }
  return(data)
}
