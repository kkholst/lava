# ----- Deprecated: removed from lava in version 1.9.1 -----
# Monte Carlo simulation-based p-values via estimate().
# This code is preserved for reference only and is NOT loaded by the package.

# ---- estimate_sim ----------------------------------------------------------
# Standalone replacement for the removed estimate(..., R=, null.sim=) pathway.
#
# @param x model object (glm, lvmfit, ...)
# @param f transformation of model parameters
# @param R number of simulations
# @param null.sim mean under the null for simulations
# @param level level of confidence limits
# @param vcov (optional) covariance matrix of parameter estimates
# @param coef (optional) parameter coefficients
# @param labels (optional) names of coefficients
# @param robust if TRUE robust standard errors are used
# @param IC if TRUE the influence function decomposition is used
# @param ... additional arguments
estimate_sim <- function(
  x = NULL,
  f,
  R,
  null.sim,
  ...,
  level = 0.95,
  vcov,
  coef,
  labels,
  robust = TRUE,
  IC = robust
) {
  if (is.null(f)) {
    stop("Supply function 'f'")
  }
  if (!missing(coef)) {
    pp <- coef
  } else {
    pp <- stats::coef(x)
  }
  if (missing(vcov) || is.null(vcov)) {
    if (!is.null(IC) && ((is.logical(IC) && IC) || length(IC) > 0) && robust) {
      ic_theta <- lava::IC(x)
      V <- var(ic_theta) / nrow(ic_theta)
    } else {
      V <- stats::vcov(x)
    }
  } else {
    V <- cbind(vcov)
  }
  if (missing(null.sim)) {
    null.sim <- rep(0, length(pp))
  }
  est <- f(pp)
  if (is.list(est)) {
    nn <- names(est)
    est <- unlist(est)
    names(est) <- nn
  }
  if (missing(labels)) {
    labels <- colnames(rbind(est))
  }
  res <- simnull(R, f, mu = null.sim, sigma = V, labels = labels)
  structure(
    res,
    class = c("estimate.sim", "sim"),
    coef = pp,
    vcov = V,
    f = f,
    estimate = est
  )
}

simnull <- function(R, f, mu, sigma, labels = NULL) {
  X <- lava::rmvn0(R, mu = mu, sigma = sigma)
  est <- f(mu)
  res <- apply(X, 1, f)
  if (is.list(est)) {
    nn <- names(est)
    est <- unlist(est)
    names(est) <- nn
    res <- matrix(unlist(res), byrow = TRUE, ncol = length(est))
  } else {
    res <- t(rbind(res))
  }
  if (is.null(labels)) {
    labels <- colnames(rbind(est))
    if (is.null(labels)) {
      labels <- paste0("p", seq_along(est))
    }
  }
  colnames(res) <- labels
  return(res)
}

estimate.estimate.sim <- function(x, f, R = 0, labels, ...) {
  atr <- attributes(x)
  if (R > 0) {
    if (missing(f)) {
      val <- simnull(
        R,
        f = atr[["f"]],
        mu = atr[["coef"]],
        sigma = atr[["vcov"]]
      )
      res <- rbind(x, val)
      for (a in setdiff(names(atr), c("dim", "dimnames"))) {
        attr(res, a) <- atr[[a]]
      }
    } else {
      res <- simnull(R, f = f, mu = atr[["coef"]], sigma = atr[["vcov"]])
      for (a in setdiff(names(atr), c("dim", "dimnames", "f"))) {
        attr(res, a) <- atr[[a]]
      }
      attr(f, "f") <- f
      est <- unlist(f(atr[["coef"]]))
      if (missing(labels)) {
        labels <- colnames(rbind(est))
      }
      attr(res, "estimate") <- est
    }
    if (!missing(labels)) {
      colnames(res) <- labels
    }
    return(res)
  }
  if (missing(f)) {
    if (!missing(labels)) {
      colnames(res) <- labels
    }
    return(x)
  }

  est <- f(atr[["coef"]])
  res <- apply(x, 1, f)
  if (is.list(est)) {
    res <- matrix(unlist(res), byrow = TRUE, ncol = length(est))
  } else {
    res <- t(rbind(res))
  }
  if (missing(labels)) {
    labels <- colnames(rbind(est))
    if (is.null(labels)) labels <- paste0("p", seq_along(est))
  }
  colnames(res) <- labels
  for (a in setdiff(names(atr), c("dim", "dimnames", "f", "estimate"))) {
    attr(res, a) <- atr[[a]]
  }
  attr(f, "f") <- f
  attr(res, "estimate") <- unlist(est)
  return(res)
}

print.estimate.sim <- function(x, level = .95, ...) {
  quantiles <- c((1 - level) / 2, 1 - (1 - level) / 2)
  est <- attr(x, "estimate")
  mysummary <- function(x, INDEX, ...) {
    x <- as.vector(x)
    res <- c(
      mean(x, na.rm = TRUE),
      sd(x, na.rm = TRUE),
      quantile(x, quantiles, na.rm = TRUE),
      est[INDEX],
      mean(abs(x) > abs(est[INDEX]), na.rm = TRUE)
    )

    names(res) <- c(
      "Mean",
      "SD",
      paste0(quantiles * 100, "%"),
      "Estimate",
      "P-value"
    )
    res
  }
  env <- new.env()
  assign("est", attr(x, "estimate"), env)
  environment(mysummary) <- env
  print(summary(x, fun = mysummary, ...))
}

# ---- examples --------------------------------------------------------------
# Monte Carlo approach, simple trend test example
#
# library(lava)
# m <- categorical(lvm(), ~x, K = 5)
# regression(m, additive = TRUE) <- y ~ x
# d <- simulate(m, 100, seed = 1, 'y~x' = 0.1)
# l <- lm(y ~ -1 + factor(x), data = d)
#
# f <- function(x) coef(lm(x ~ seq_along(x)))[2]
# null <- rep(mean(coef(l)), length(coef(l)))
# ## just need to make sure we simulate under H0: slope = 0
# source(system.file("misc/estimate_sim.R", package = "lava"))
# estimate_sim(l, f, R = 1e2, null.sim = null)
