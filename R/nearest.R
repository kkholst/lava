##' Find nearest element
##' 
##' Find elements in a vector that are close to a numeric \code{x}
##' @param x \code{numeric}
##' @param yy vector of \code{numeric}
##' @param epsilon threshold
##' @param norm Function that measures the distance
##' @return Returns the indices of which elements in \code{yy} that satisfy
##' \deqn{||y-x||<\epsilon}.
##' @author Klaus K. Holst
##' @keywords utilities
##' @export
`nearest` <-
function(x,yy,epsilon,norm=function(z) z^2) { 
  if (!is.numeric(x) & !is.numeric(yy)) stop("numeric variables only!")
  dif <- norm(yy-x)
  if (is.missing(epsilon)) return(which.min(dif)) 
  return(which(dif<epsilon))
}


