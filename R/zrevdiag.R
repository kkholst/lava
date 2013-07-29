##' @title Create/extract 'reverse'-diagonal matrix
##' @aliases revdiag revdiag<-
##' @usage
##' revdiag(x,...)
##' 
##' revdiag(x,...) <- value
##' @param x vector
##' @param value For the assignment function the values to put in the diagonal
##' @param \dots additional arguments to lower level functions
##' @author Klaus K. Holst
##' @export
revdiag <- function(x,...) {
    if (NCOL(x)==1) {
      res <- matrix(0,length(x),length(x))
      revdiag(res) <- x
      return(res)
    }
    n <- ncol(x)
    x[cbind(rev(seq(n)),seq(n))]
  }

##' @export
"revdiag<-" <- function(x,...,value) {
  n <- ncol(x)
  x[cbind(rev(seq(n)),seq(n))] <- value
  x
}
