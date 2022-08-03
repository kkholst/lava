##' Extract i.i.d. decomposition from model object
##'
##' This function extracts
##' @param x Model object
##' @param ... Additional arguments (see the man-page of the IC method)
##' @export
iid <- function(x, ...) UseMethod("iid")

##' @export
iid.default <- function(x, ...) {
  res <- IC(x,...)
  if (!is.null(attr(res, "bread"))) {
    attr(res, "bread") <- attr(res, "bread")/NROW(res)
  }
  return(res/NROW(res))
}

