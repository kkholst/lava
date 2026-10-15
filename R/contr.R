#' Create contrast matrix
#'
#' Create contrast matrix typically for use with 'estimate' (Wald tests).
#' @export
#' @param p index of non-zero entries (see example)
#' @param n Total number of parameters (if omitted the max number in p will be
#'   used)
#' @param diff If FALSE all non-zero entries are +1, otherwise the second
#'   non-zero element in each row will be -1.
#' @param ... Additional arguments to lower level functions
#' @aliases contr parsedesign
#' @seealso [pairwise_diff]
#' @examples
#' contr(2,n=5)
#' contr(as.list(2:4),n=5)
#' contr(list(1,2,4),n=5)
#' contr(c(2,3,4),n=5)
#' contr(list(c(1,3),c(2,4)),n=5)
#' contr(list(c(1,3),c(2,4),5))
#'
#' parsedesign(c("aa","b","c"),"?","?",diff=c(FALSE,TRUE))
#'
#' ## All pairs comparisons:
#' pdiff <- function(n) lava::contr(lapply(seq(n-1), function(x) seq(x, n)))
#' pdiff(4)
contr <- function(p, n, diff = TRUE, ...) {
  if (missing(n) || is.null(n)) {
    n <- max(unlist(p))
  }
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

#' @title Pairwise parameter differences
#'
#' @description
#' Construct contrasts for pairwise differences between parameters,
#' or apply these contrasts to an `estimate` object.
#'
#' @param n An integer giving the total number of parameters, or an
#'   `estimate` object from which the number of parameters is derived.
#' @param idx An integer vector of parameter indices to compare.
#'   If `NULL`, all parameters are compared.
#' @param reverse Logical; if `TRUE`, reverse the direction of each
#'   difference by multiplying the contrast matrix by -1.
#'
#' @details
#' When `n` is an integer, the function returns a contrast matrix with
#'   `n` columns and one row for each pairwise comparison. If `idx` is
#'   supplied, only pairs of parameters selected by `idx` are compared.
#'
#' When `n` is an `estimate` object, the number of parameters is obtained
#'   from `length(coef(n))`. The contrast matrix is applied to the object,
#'   returning an `estimate` object containing the pairwise differences
#'   and their propagated uncertainty.
#'
#' @return
#' A numeric contrast matrix if `n` is an integer, or an `estimate`
#'   object if `n` is an `estimate` object.
#'
#' @examples
#' pairwise_diff(4)
#' pairwise_diff(4, idx = c(1, 3, 4))
#' pairwise_diff(4, reverse = TRUE)
#'
#' e <- estimate(coef = c(a = 1, b = 2, c = 4),
#'               vcov = diag(3))
#' pairwise_diff(e)
#' @export
#' @aliases pairwise_diff pairwise.diff
pairwise_diff <- function(n, idx=NULL, reverse=FALSE) {
  estimate_obj <- inherits(n, "estimate")
  if (isTRUE(estimate_obj)) {
    e <- n
    n <- length(coef(e))
  }
  if (!is.null(idx)) {
    if (any(idx > n)) {
      stop("'n' must be larger or equal to any element of 'idx'")
    }
    B <- contr(combn(idx, 2, simplify=FALSE), n)
  } else {
    pdiff <- function(n) lava::contr(lapply(seq(n-1), function(x) seq(x, n)))
    B <- pdiff(n)
  }
  if (reverse) B <- -B
  if (isTRUE(estimate_obj)) {
    return(B %*% e)
  }
  return(B)
}

##' @export
pairwise.diff <- function(n) {
  .Deprecated("pairwise_diff", package = "lava",
              msg = paste(
                "`pairwise.diff` is deprecated. Use `pairwise_diff` instead."
              )
              )
  pairwise_diff(n)
}
