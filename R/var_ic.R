#' @export
var_ic <- function(x, ...) {
  N <- crossprod(!is.na(x))
  V <- var(x, use="pairwise.complete.obs")*(N-1)/N^2
  V[N==0] <- 0
  return(V)
}
