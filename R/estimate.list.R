##' Estimate method for lists
##'
##' @param x list of model objects or `lvm` objects
##' @param ... Additional arguments to lower level functions
##' @return Object of class `estimate.list` (a list of `estimate` objects),
##'   or a fitted multi-group model if input is a list of `lvm` objects.
#' @export
estimate.list <- function(x, ...) {
  if (inherits(x[[1]], "lvm")) {
    return(estimate_lvmlist(x, ...))
  }
  res <- lapply(x, function(x) estimate(x, ...))
  class(res) <- c("estimate.list", "list")
  res
}

#' @export
coef.estimate.list <- function(object, ...) {
  lapply(object, coef)
}

#' @export
vcov.estimate.list <- function(object, ...) {
  lapply(object, vcov)
}

#' @export
IC.estimate.list <- function(x, ...) {
  lapply(x, vcov)
}

#' @export
merge.estimate.list <- function(x, ...) {
  Reduce(merge, x)
}
