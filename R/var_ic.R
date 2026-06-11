##' Variance based on influence function
##'
##' Computes the empirical variance of the influence function, providing a
##' sandwich-type variance estimator.
##' @param x Influence function matrix (observations x parameters), typically
##'   obtained from [IC].
##' @param ... Additional arguments (currently not used)
##' @return Covariance matrix (p x p).
#' @export
var_ic <- function(x, ...) {
  N <- crossprod(!is.na(x))
  V <- var(x, use="pairwise.complete.obs")*(N-1)/N^2
  V[N==0] <- 0
  return(V)
}
