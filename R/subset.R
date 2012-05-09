##' Extract subset of latent variable model
##' 
##' Extract measurement models or user-specified subset of model
##' 
##' 
##' @aliases measurement
##' @param x \code{lvm}-object.
##' @param vars Character vector or formula specifying variables to include in
##' subset.
##' @param \dots Additional arguments to be passed to the low level functions
##' @return A \code{lvm}-object.
##' @author Klaus K. Holst
##' @keywords models regression
##' @examples
##' 
##' m <- lvm(c(y1,y2)~x1+x2)
##' subset(m,~y1+x1)
##'
##' @S3method subset lvm
##' @method subset lvm
subset.lvm <- function(x, vars, ...) {
  if (missing(vars))
    return(x)
  if (class(vars)[1]=="formula") vars <- all.vars(vars)
  if (!all(vars%in%vars(x))) stop("Not a subset of model")  
  latentvars <- intersect(vars,latent(x))
  g0 <- subGraph(vars, Graph(x))
  res <- graph2lvm(g0)
  if (length(latentvars)>0)
    latent(res) <- latentvars
  res$cov[vars,vars] <- x$cov[vars,vars]
  ## Fixed parameters:
  res$par[vars,vars] <- x$par[vars,vars]
  res$fix[vars,vars] <- x$fix[vars,vars]
  res$covpar[vars,vars] <- x$covpar[vars,vars]
  res$covfix[vars,vars] <- x$covfix[vars,vars]
  res$mean[vars] <- x$mean[vars]
  index(res) <- reindex(res)
  return(res)
}
