##' @export
`Weight` <- function(x,...) UseMethod("Weight")

##' @S3method Weight default
Weight.default <- function(x,...) eval(x$weight)
