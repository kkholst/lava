"kill<-" <- function(x, ..., value) UseMethod("kill<-")

"kill<-.lvm" <- function(x, ..., value) {
  if (class(value)[1]=="formula") {
    return(kill(x,all.vars(value)))
  }  
  kill(x,value)
}

"kill" <- function(x, value, ...) UseMethod("kill")
"kill.lvm" <- function(x, value, ...) {
  idx <- which(vars(x)%in%value)
  if (length(idx)==0)
    return(x)
  vv <- vars(x)[idx]
  keep <- setdiff(1:length(vars(x)),idx)
  M <- as(Graph(x), Class="matrix")
  for (v1 in vv)
    Graph(x) <- removeNode(v1, Graph(x))
  x$par <- x$par[keep,keep,drop=FALSE]
  x$fix <- x$fix[keep,keep,drop=FALSE]
  x$covpar <- x$covpar[keep,keep,drop=FALSE]
  x$covfix <- x$covfix[keep,keep,drop=FALSE]
  x$cov <- x$cov[keep,keep,drop=FALSE]
  x$mean <- (x$mean)[-idx]
  x$exogenous <- setdiff(exogenous(x),vv)
  index(x) <- reindex(x)
  return(x)
}
