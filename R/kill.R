#' Add or remove variables from latent variables. Remove associations from
#' latent variable model
#' 
#' Add or remove variables from a latent variable model, or remove previously
#' specified associations among variables.
#' 
#' 
#' @aliases kill kill<- kill.lvm cancel cancel cancel<- cancel.lvm addvar<-
#' @param x \code{lvm}-object
#' @param value Vector of variables or formula specifying which nodes to
#' remove/add (or associations between nodes).
#' @param \dots ...
#' @return A \code{lvm}-object
#' @author Klaus K. Holst
#' @keywords models regression
#' @examples
#' 
#' m <- lvm()
#' addvar(m) <- ~y1+y2+x
#' covariance(m) <- y1~y2
#' regression(m) <- c(y1,y2) ~ x
#' ### Cancel the covariance between the residuals of y1 and y2
#' cancel(m) <- y1~y2
#' ### Remove y2 from the model 
#' kill(m) <- ~y2
#' 
"kill<-" <- function(x, ..., value) UseMethod("kill<-")

"kill<-.lvm" <- function(x, ..., value) {
  kill(x,value)
}

"kill" <- function(x, value, ...) UseMethod("kill")
"kill.lvm" <- function(x, value, ...) {
  if (class(value)[1]=="formula") {
    return(kill(x,all.vars(value),...))
  }  
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
