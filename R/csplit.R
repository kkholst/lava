##' Split data into folds
##'
##' @title Split data into folds
##' @param x Data or integer (size)
##' @param p Number of folds, or if a number between 0 and 1 is given two folds of size p and (1-p) will be returned
##' @param replace With or with-out replacement
##' @param return.index If TRUE index of folds are returned otherwise the actual data splits are returned (default)
##' @param k (Optional, only used when p=NULL) number of folds without shuffling
##' @param ... additional arguments to lower-level functions
##' @export
##' @aliases csplit foldr
##' @examples
##' foldr(5,2,rep=2)
##' csplit(10,3)
##' csplit(iris[1:10,]) ## Split in two sets 1:(n/2) and (n/2+1):n
##' csplit(iris[1:10,],0.5)
##' @author Klaus K. Holst
csplit <- function(x, p=NULL, replace=FALSE,
                   return.index=FALSE, k=2, ...) {
  if (length(x)==1 && is.numeric(x))
    x <- seq(x)
  N <- NROW(x)
  if (is.null(p)) { ##
    K <- base::round(N/k)
    idx <- split(seq(N), sort(rep(seq(k), length.out=N, each=K)))
  } else {
    if (p<1) { ## two folds (size N*p and N*(1-p))

      idx1 <- base::sample(N, base::round(p*N), replace=replace)
      idx <- list(idx1,
                  base::sample(setdiff(seq(N), idx1), replace=replace))
    } else { ## Number of folds (equal size)
      idx <- split(sample(seq(N)), rep(seq(p), length=N))
    }
  }
  if (return.index)
    return(idx)
  if (!is.vector(x)) {
    return(lapply(idx, function(ii) x[ii, , drop=FALSE]))
  }
  return(lapply(idx, function(ii) x[ii]))
}

##' @export
foldr <- function(n, K=5, rep=1, list=TRUE) {
  res <- replicate(rep, split(sample(seq(n)),
                              rep(seq(K), length=n)),
                   simplify=FALSE)
  if (!list) {
    ids <- rep(seq_len(length(res[[1]])),
           unlist(lapply(res[[1]], length)))
    res <- lapply(res, function(x) {
      ids[order(unlist(x))]
    })
    if (length(res)==1)
      res <- res[[1]]
  }
  return(res)
}

