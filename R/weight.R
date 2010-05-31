`Weight` <- function(x,...) UseMethod("Weight")
Weight.default <- function(x,...) eval(x$weight)
