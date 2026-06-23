with_unique_warnings <- function(expr) {
  # utility function that prevents the same warning being cast more than once
  seen <- character(0)
  withCallingHandlers(expr, warning = function(w) {
    msg <- conditionMessage(w)
    if (msg %in% seen) {
      invokeRestart("muffleWarning")
    } else {
      seen <<- c(seen, msg)
    }
  })
}


small_sample_correction <- function(ic, type, var.adj) {
  V <- var_ic(ic)
  ## Small-sample corrections for clustered data
  K <- NROW(ic)
  N <- attributes(ic)$N
  if (is.null(N)) N <- K
  p <- NCOL(ic)
  adj0 <- K/(K-p) ## Mancl & DeRouen, 2001
  adj1 <- K/(K-1) ## Mancl & DeRouen, 2001
  adj2 <- (N-1)/(N-p)*(K/(K-1)) ## Morel,Bokossa & Neerchal, 2003
  if (tolower(type[1])=="mbn" && !is.null(attributes(ic)$bread)) {
    V0 <- V
    iI0 <- attributes(ic)$bread
    I0 <- Inverse(iI0)
    delta <- min(0.5, p / (K - p))
    phi <- max(1, tr(I0%*%V0)*adj2/p)
    V <- adj2*V0 + delta*phi*iI0
  }
  if (tolower(type[1])=="df") {
    V <- adj0*V
  }
  if (tolower(type[1])=="df1") {
    V <- adj1*V
  }
  if (tolower(type[1])=="df2") {
    V <- adj2*V
  }
  if (tolower(type[1])%in%c("hc3", "hc4")) {
    ic <- cbind(ic)
    S <- Inverse(crossprod(ic), tol=sqrt(.Machine$double.eps))
    h_emp <- rowSums((ic %*% S) * ic) # empirical h, lev.
    n <- nrow(ic)
    ## h_emp <- pmin(h_emp, 0.99) * (n-1)/n  # Truncate leverage to prevent division by zero
    if (tolower(type[1])=="hc3") {
      ## phi_norm <- sqrt(rowSums(ic^2))
      ## ex_kurt <- (mean((phi_norm - mean(phi_norm))^4) / var(phi_norm)^2) - 3
      ## alpha <- exp(-max(0, ex_kurt) / 25)
      v.alpha <- ifelse(missing(var.adj), 0.25, var.adj)
      h <- v.alpha * h_emp + (1 - v.alpha) * (ncol(ic) / n)
      adj <- 1 / (1 - h)
    } else {
      ## Cribari-Neto (2004)
      delta <- pmin(1.5, n/ncol(ic) * h_emp)
      adj <- 1 / (1 - h_emp)**delta
    }
    for (i in seq_len(NCOL(ic))) ic[, i] <- ic[, i] * adj
    V <- var_ic(ic)
  }
  return(V)
}

# construct coefficient matrix with confidence limits and two-sided p-values
estimate_coefmat <- function(est, se, df, level = 0.95, null = 0) {
  alpha <- 1 - level
  alpha.str <- paste(c(alpha/2, 1 -alpha/2)*100, "", sep="%")
  if (!is.null(df)) {
    za <- qt(1-alpha/2, df=df)
    pval <- 2*pt(abs((est-null) / se), df=df, lower.tail=FALSE)
  } else {
    za <- qnorm(1-alpha/2)
    pval <- 2*pnorm(abs((est-null)/ se), lower.tail=FALSE)
  }
  res <- cbind(est, se, est - za * se, est + za * se, pval)
  colnames(res) <- c("Estimate", "Std.Err", alpha.str, "P-value")
  return(res)
}

#' Summary of estimate objects
#'
#' Computes hypothesis tests, contrasts, and small-sample corrections for
#' an [estimate] object. The arguments `null`, `contrast`, `type`, and
#' `var.adj` were previously available on [estimate.default()] and have
#' been moved here.
#'
#' @param object an `estimate` object.
#' @param contrast (optional) contrast matrix for the final Wald test.
#' @param null (optional) null hypothesis to test.
#' @param type type of small-sample correction. Requires the estimate
#'   to have been computed with `IC=TRUE` (the default).
#' @param var.adj variance adjustment parameter for small-sample
#'   correction. Requires the estimate to have been computed with
#'   `IC=TRUE` (the default).
#' @param df degrees of freedom for t-based inference (default: `NULL` for
#'   Gaussian approximation; when set, confidence intervals and p-values use
#'   the t-distribution with `df` degrees of freedom)
#' @param level level of confidence limits (default 0.95)
#' @param transform (optional) function applied to the point estimates
#'   and confidence interval bounds *after* inference is performed on the
#'   original scale. Useful for variance-stabilizing transformations, e.g.,
#'   compute CIs on the `atanh` (Fisher z) scale and back-transform with
#'   `tanh`.
#' @param ... additional arguments passed to [estimate()].
#' @seealso [estimate.default()]
#' @export
summary.estimate <- function(object,
                             contrast,
                             null = 0,
                             level = 0.95,
                             type,
                             var.adj = 0.25,
                             df,
                             transform,
                             ...) {
  with_unique_warnings({

    p <- coef(object)
    if (missing(df)) df <- object$df
    if (missing(contrast)) contrast <- diag(1, nrow=length(p))
    correction <- !missing(type)
    if (correction) {
      if (is.null(object$IC)) {
        stop(
          "Small-sample corrections in summary() require an estimate ",
          "computed with IC=TRUE."
        )
      }
      V <- small_sample_correction(IC(object), type=type, var.adj=var.adj)
    } else {
      V <- vcov(object)
    }
## Preserve dimnames on vcov: the recall path through type/var.adj
  ## corrections builds V from the influence function and may drop
  ## dimnames. Restore them from the original coef names so the
  ## summary output is consistent with the deprecated estimate() path.
  if (!is.null(V) && is.null(dimnames(V)) &&
      !is.null(names(p)) &&
      nrow(V) == length(p)) {
    dimnames(V) <- list(names(p), names(p))
  }

  if (is.vector(contrast) || is.list(contrast)) {
    contrast <- contr(contrast, names(res$coef), ...)
  }

  cc0 <- estimate_coefmat(p, diag(V)**.5, df=df, level=level, null=null)
  rownames(cc0) <- rownames(parameter(object))
  waldtest <- compare(object, contrast=contrast, null=null, vcov=V)
  class(object) <- "list"
  res <- c(object[c("coef", "coefmat", "vcov", "call",
                    "ncluster", "model.index")], list(compare=waldtest))
  res$coefmat <- cc0
  if (!missing(transform)) {
    res$coefmat[, c(1, 3, 4)] <- do.call(transform,
                                         list(res$coefmat[, c(1, 3, 4)]))
    res$coefmat[, 2] <- NA
  }
  class(res) <- "summary.estimate"
  return(res)
  })
}


#' @export
print.summary.estimate <- function(x, ...) {
  print.estimate(x, type=2L, ...)
}

#' @export
coef.summary.estimate <- function(object, ...) {
  return(object$coef)
}

#' @export
parameter.summary.estimate <- function(x, ...) {
  return(x$coefmat)
}

#' @export
vcov.summary.estimate <- function(object, ...) {
  return(object$vcov)
}
