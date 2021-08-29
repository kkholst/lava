##' @export
"intervention<-" <- function(object, ..., value)
  UseMethod("intervention<-")

##' @export
`intervention` <-
  function(object, ...) UseMethod("intervention")

##' Define intervention
##'
##' Define intervention in a `lvm` object
##' @param object lvm object
##' @param to String defining variable or formula
##' @param value function defining intervention
##' @param dist Distribution
##' @param ... Additional arguments to lower level functions
##' @aliases intervention<- intervention intervention.lvm intervention<-.lvm
##' @seealso regression lvm sim
##' @examples
##' m <- lvm(y ~ a + x, a ~ x)
##' distribution(m, ~a+y) <- binomial.lvm()
##' mm <- intervention(m, "a", value=3)
##' sim(mm, 10)
##' mm <- intervention(m, a~x, function(x) (x>0)*1)
##' sim(mm, 10)
##' @export
intervention.lvm <- function(object, to, value, dist=none.lvm(), ...) {
  if (!is.numeric(value))
    regression(object, to, ...) <- value
  y <- to
  if (inherits(to, "formula")) {
    y <- getoutcome(to)
    if (length(y)==0)
      y <- attr(y, "x")
  }
  parents <- parents(object, y)
  if (length(parents)>0)
    cancel(object) <- toformula(y, parents)
  if (is.numeric(value)) {
    distribution(object, y) <- constant.lvm(value)
  } else {
    distribution(object, y) <- dist
  }
  return(object)
}

##' @export
"intervention<-.lvm" <- function(object, to, ..., value) {
  object <- intervention(object, to, value, ...)
  return(object)
}
