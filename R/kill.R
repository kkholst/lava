##' Generic method for removing elements of object
##'
##' @title Remove variables from (model) object. 
##' @aliases kill kill<-
##' @param x Model object
##' @param value Vector of variables or formula specifying which nodes to
##' remove/add (or associations between nodes).
##' @param \dots ..
##' @usage
##' kill(x, ...) <- value
##' @seealso \code{cancel}
##' @author Klaus K. Holst
##' @keywords models regression
##' @export
##' @examples
##' 
##' m <- lvm()
##' addvar(m) <- ~y1+y2+x
##' covariance(m) <- y1~y2
##' regression(m) <- c(y1,y2) ~ x
##' ### Cancel the covariance between the residuals of y1 and y2
##' cancel(m) <- y1~y2
##' ### Remove y2 from the model 
##' kill(m) <- ~y2
##'
"kill" <- function(x, ...) UseMethod("kill")

##' @export
"kill<-" <- function(x, ..., value) UseMethod("kill<-")

##' @S3method kill<- lvm
"kill<-.lvm" <- function(x, ..., value) {
  kill(x,value)
}

##' @S3method kill lvm
"kill.lvm" <- function(x, value, ...) {
  if (class(value)[1]=="formula") {
    return(kill(x,all.vars(value),...))
  }  
  idx <- which(vars(x)%in%value)
  if (length(idx)==0)
    return(x)
  vv <- vars(x)[idx]
  keep <- setdiff(1:length(vars(x)),idx)
  x$M <- x$M[keep,keep,drop=FALSE]
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
