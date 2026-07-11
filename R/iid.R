##' Extract i.i.d. decomposition from model object
##'
##' This function extracts
##' @param x Model object
##' @param ... Additional arguments (see the man-page of the IC method)
##' @return Matrix (n x p) with the i.i.d. decomposition of the estimator
##'   (scaled influence function, i.e., the influence function divided by n).
##' @export
iid <- function(x, ...) UseMethod("iid")

##' @export
iid.default <- function(x, ...) {
  res <- IC(x, ...)
  if (!is.null(attr(res, "bread"))) {
    attr(res, "bread") <- attr(res, "bread") / NROW(res)
  }
  return(res / NROW(res))
}
