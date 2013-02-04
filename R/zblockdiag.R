##' @title Combine matrices to block diagonal structure
##' @param x Matrix
##' @param \dots Additional matrices
##' @param pad Vyalue outside block-diagonal
##' @author Klaus K. Holst
##' @export
##' @examples
##' A <- diag(3)+1
##' blockdiag(A,A,A,pad=NA)
blockdiag <- function(x,...,pad=0) {
  if (is.list(x)) xx <- x  else xx <- list(x,...)
  xx <- list(x,...)
  rows <- unlist(lapply(xx,nrow))
  crows <- c(0,cumsum(rows))
  cols <- unlist(lapply(xx,ncol))
  ccols <- c(0,cumsum(cols))
  res <- matrix(pad,nrow=sum(rows),ncol=sum(cols))
  for (i in 1:length(xx)) {
    idx1 <- 1:rows[i]+crows[i]; idx2 <- 1:cols[i]+ccols[i]
    res[idx1,idx2] <- xx[[i]]
  }
  colnames(res) <- unlist(lapply(xx,colnames)); rownames(res) <- unlist(lapply(xx,rownames))
  return(res)
}

