##' Finds elements in vector or column-names in data.frame/matrix
##' 
##' Pattern matching in a vector or column names of a data.frame or matrix.
##' @param x vector, matrix or data.frame.
##' @param pattern regular expression to search for
##' @param subset If TRUE returns subset of data.frame/matrix otherwise just the matching column names
##' @param ignore.case Default ignore case
##' @param ... 
##' @return A data.frame with 2 columns with the indices in the first and the
##' matching names in the second.
##' @author Klaus K. Holst
##' @seealso \code{\link{grep}}, and \code{\link{agrep}} for approximate string
##' matching,
##' @keywords misc utilities
##' @examples
##' data(iris)
##' head(Grep(iris,"(len)|(sp)"))
##' @export
`Grep` <-
function(x, pattern, subset=TRUE, ignore.case = TRUE,...) {
  if (is.data.frame(x))
    nn <- names(x)
  else if (is.matrix(x))
      nn <- colnames(nn)
  else nn <-  x
  ii <- grep(pattern,nn,ignore.case=ignore.case,...)
  if (subset) {
      if (is.matrix(x) || is.data.frame(x))
          return(x[,ii,drop=FALSE])
      else return(x[ii])
  }  
  res <- data.frame(index=ii,name=nn[ii]);
  res
}

