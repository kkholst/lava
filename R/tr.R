##' Trace operator
##'
##' Calculates the trace of a square matrix.
##' @param x Square numeric matrix
##' @param \dots Additional arguments to lower level functions
##' @return \code{numeric}
##' @author Klaus K. Holst
##' @seealso \code{\link{crossprod}}, \code{\link{tcrossprod}}
##' @keywords math algebra
##' @examples
##'
##' tr(diag(1:5))
##' @export
"tr" <- function(x,...) UseMethod("tr")

##' @export
`tr.matrix` <-
function(x,na.rm=FALSE,...) {
  if (length(x)==1)
    return(x)
  n <- nrow(x)
  if (!n)
    stop("0 x 0 matrix")
  if (n != ncol(x))
    stop("non-square matrix")
  if (!na.rm && any(!is.finite(x)))
    stop("infinite or missing values")
  return(sum(diag(x),na.rm=na.rm))
}
