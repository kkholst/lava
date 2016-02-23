##' Create/extract 'reverse'-diagonal matrix
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

##' @export
offdiag <- function(x,...) {
    if (NCOL(x)==1) return(NULL)
    ii <- c(which(lower.tri(x,diag=FALSE)),which(upper.tri(x,diag=FALSE)))
    res <- x[ii]
    class(res) <- c("offdiag",class(res))
    res
  }

##' @export
"offdiag<-" <- function(x,...,value) {
    ii <- c(which(lower.tri(x,diag=FALSE)),which(upper.tri(x,diag=FALSE)))   
    x[ii] <- value
    return(x)
}

##' @export
print.offdiag <- function(x,...) {
    n <- (1+sqrt(1+4*length(x)))/2
    M <- matrix(NA,n,n)
    M[lower.tri(M)] <- x[seq(length(x)/2)]
    M[upper.tri(M)] <- x[seq(length(x)/2)+length(x)/2]
    print(M,na.print="")
}
