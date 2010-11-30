
`heavytail` <- function(x,...) UseMethod("heavytail")
"heavytail<-" <- function(x,...,value) UseMethod("heavytail<-")

"heavytail<-.lvm" <- function(x,...,value) {
  if (class(value)[1]=="formula") {
    return(heavytail(x,all.vars(value),...))
  }
  heavytail(x, value, ...)
}

`heavytail.lvm` <-
function(x,var=NULL,df=1,...) {
  if (is.null(var)) {
    htidx <- unlist(nodeData(Graph(x), attr="heavytail"))
    if (length(htidx)>0 && any(htidx!=0)) {
      res <- htidx[htidx>0]
      attributes(res)$couple <- unlist(nodeData(Graph(x), attr="heavytail.couple"))[htidx>0]
      return(res)      
    }
    return(NULL)
  }
  couples <- attributes(heavytail(x))$couple
  if (length(couples)==0) newval <- 1
  else newval <- max(couples)+1
  nodeData(Graph(x), var, attr="heavytail.couple") <- newval
  nodeData(Graph(x), var, attr="heavytail") <- df
  return(x)
}

heavytail.init.hook <- function(x,...) {
  nodeDataDefaults(x,"heavytail") <- 0
  nodeDataDefaults(x,"heavytail.couple") <- 0
  return(x)
}
heavytail.sim.hook <- function(x,data,...) {
  n <- nrow(data)
  hvar <- heavytail(x)
  if (length(hvar)==0) return(data)
  couples <- unique(attributes(hvar)$couple)
  h.type <- list()
  for (j in couples)
    h.type <- c(h.type, list( hvar[(which(attributes(hvar)$couple==j))]))
  for (i in 1:length(couples)) {
    df <- hvar[[i]][1]
    Z <- rchisq(n,df=df)/df
    for (v in names(h.type[[1]])) 
      data[,v] <- data[,v]/sqrt(Z)
  }
  return(data)
}
