#' @export
estimate.data.frame <- function(x, ...) {
  estimate(as.matrix(x), ...)
}

IC_quantile <- function(x, estimate, probs=0.5, ...) {
  x <- na.omit(x)
  f0 <- density(x, ...)
  ## U <- function(est) (tau - (x <= est))
  if (missing(estimate)) {
    estimate <- quantile(x, probs=probs)
  }
  res <- c()
  for (i in seq_len(length(estimate))) {
    res <- cbind(res,
    (probs[i] - (x <= estimate[i]))/with(f0, approx(x, y, estimate[i]))$y
    )
  }
  res
}

#' Estimate parameters and influence function.
#'
#' Estimate parameters for the sample mean, variance, and quantiles
#' @export
#' @aliases estimate.array estimate.data.frame
#' @param x numeric matrix
#' @param type target parameter ("mean", "variance", "quantile")
#' @param probs numeric vector of probabilities (for type="quantile")
#' @param ... Additional arguments to lower level functions (i.e.,
#'   stats::density.default when type="quantile")
#' @return Object of class `estimate` (see [estimate.default]).
estimate.array <- function(x, type="mean", probs=0.5, ...) {
  cl <- match.call()
  if (missing(x) || is.null(x)) {
    return(estimate(NULL, ...))
  }
  dots <- list(...)
  density.args <- dots[]
  cc <- apply(x, 2, function(y) mean(y, na.rm = TRUE))
  ic <- apply(x, 2, function(y) y - mean(y, na.rm = TRUE))
  if (tolower(type) %in% c("var", "variance")) {
    cc <- apply(x, 2, function(y) mean((y - mean(y, na.rm=TRUE))^2, na.rm = TRUE))
    ic <- ic^2
    for (i in seq_len(NCOL(ic))) {
      ic[, i] <- ic[, i] - cc[i]
    }
  }
  if (tolower(type) %in% c("quantile")) {
    density.args <- list()
    dargs <- names(formals(density.default))
    didx <- which(dargs %in% names(dots))
    if (length(didx)>0) {
      density.args <- dots[dargs[didx]]
      dots[dargs[didx]] <- NULL
    }
    cc <- unlist(apply(x, 2, function(y)
      quantile(y, probs=probs, na.rm = TRUE),
      simplify=FALSE))
    ic <- c()
    for (i in seq_len(NCOL(x))) {
      ic <- cbind(ic, do.call(IC_quantile,
                        c(list(x[, i], probs=probs), density.args)))
    }
  }
  if (any(c("vcov", "IC") %in% names(list(...)))) {
    return(estimate(NULL, coef = cc, ...))
  }
  res <- do.call(estimate, c(list(NULL, coef = cc, IC = ic), dots))
  res$call <- cl
  return(res)
}
