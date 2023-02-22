##' Create contrast matrix
##'
##' Create contrast matrix typically for use with 'estimate' (Wald tests).
##' @export
##' @param p index of non-zero entries (see example)
##' @param n Total number of parameters (if omitted the max number in p will be used)
##' @param diff If FALSE all non-zero entries are +1, otherwise the second non-zero element in each row will be -1.
##' @param ... Additional arguments to lower level functions
##' @aliases contr parsedesign pairwise.diff
##' @examples
##' contr(2,n=5)
##' contr(as.list(2:4),n=5)
##' contr(list(1,2,4),n=5)
##' contr(c(2,3,4),n=5)
##' contr(list(c(1,3),c(2,4)),n=5)
##' contr(list(c(1,3),c(2,4),5))
##'
##' parsedesign(c("aa","b","c"),"?","?",diff=c(FALSE,TRUE))
##'
##' ## All pairs comparisons:
##' pdiff <- function(n) lava::contr(lapply(seq(n-1), \(x) seq(x, n)))
##' pdiff(4)
contr <- function(p, n, diff = TRUE, ...) {
  if (missing(n)) n <- max(unlist(p))
  if (is.character(p)) {
    return(parsedesign(n, p, ...))
  }
  if (is.list(p)) {
    return(Reduce(rbind, lapply(p, function(x) {
      do.call(contr, list(x, n, diff[1L]))
    })))
  }
  if (is.character(n)) n <- length(n)
  if (!is.numeric(n)) {
    try(n <- length(coef(n)), silent = TRUE)
  }
  B <- matrix(0, ncol = n, nrow = max(1L, length(p) - 1L))
  B[, p[1]] <- 1L
  if (length(p) > 1L) {
    B[cbind(seq_len(nrow(B)), p[-1])] <- ifelse(diff[1L], -1, 1)
  }
  return(B)
}

##' @export
pairwise.diff <- function(n) {
  pdiff <- function(n) lava::contr(lapply(seq(n-1), \(x) seq(x, n)))
  pdiff(n)
}
