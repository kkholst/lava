##' @export
"intervention<-" <- function(object, ..., value)
  UseMethod("intervention<-")

##' @export
`intervention` <-
  function(object, ...) UseMethod("intervention")


##' @export
intervention.lvm <- function(object, to, value, dist=none.lvm(), ...) {
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
  distribution(object, y) <- dist
  return(object)
}

##' @export
"intervention<-.lvm" <- function(object, to, ..., value) {
  object <- intervention(object, to, value, ...)
  return(object)
}
