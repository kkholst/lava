### Return position of variance elements in the parameter vector (without mean parameters)
### Optimization constraints are needed on these parameters
variances <- function(x,mean=FALSE) {
##  if (is.null(x$parpos))
##    x$parpos <- parpos(x)
  x$parpos <- parpos(Model(x),mean=TRUE) 
  res <- diag(x$parpos$P)[diag(index(x)$P0)==1]
  if (!mean) {
    res - index(x)$npar.mean
  }
}
## And the off-diagonal (covariance) parameters
offdiags <- function(x,mean=FALSE) {
  parpos <- parpos(x,mean=mean)
  pp <- parpos$P
  pp[lower.tri(pp)][(index(x)$P0)[lower.tri(pp)]==1]
}

