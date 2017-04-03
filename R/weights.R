##' @export
`Weights` <- function(x,...) UseMethod("Weights")

##' @export
Weights.default <- function(x,...) eval(x$weights)
