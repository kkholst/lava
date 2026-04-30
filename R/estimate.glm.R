##' @export
estimate.glm <- function(x, ...) {
  estimate.default(x, ...)
}

##' @export
IC.mlm <- function(x, ...) {
  cc <- coef(x)
  r <- residuals(x)
  X <- model.matrix(x)
  w <- weights(x)
  if (!is.null(w)) {
    r <- apply(r, 2, function(x) x * w)
  }
  q <- NCOL(cc)
  ics <- lapply(1:q, function(i) {
    apply(X, 2, function(x) x * r[, i])
  })
  res <- Reduce("cbind", ics)
  colnames(res) <- names(pars(x))
  return(res)
}

##' @export
pars.mlm <- function(x, ...) {
  cc <- coef(x)
  q <- NCOL(cc)
  nn <- unlist(lapply(
    1:q, function(i)
      paste0(colnames(cc)[i], ":", rownames(cc))
  ))
  coefs <- unlist(lapply(1:q, function(x) cc[, x, drop=TRUE]))
  names(coefs) <- nn
  coefs
}

##' @export
estimate.mlm <- function(x, ...) {
  estimate.default(x, coef=pars(x), ...)
}
