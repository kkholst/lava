##' Generic method for calculating bootstrap statistics
##'
##' @title Generic bootstrap method
##' @param x Model object
##' @param \dots Additional arguments
##' @seealso \code{bootstrap.lvm} \code{bootstrap.lvmfit}
##' @author Klaus K. Holst
##' @export
bootstrap <- function(x,...) UseMethod("bootstrap")

##' Calculate bootstrap estimates of a lvm object
##'
##' Draws non-parametric bootstrap samples
##'
##' @param x \code{lvm}-object.
##' @param R Number of bootstrap samples
##' @param data The data to resample from
##' @param fun Optional function of the (bootstrapped) model-fit defining the
##' statistic of interest
##' @param control Options to the optimization routine
##' @param p Parameter vector of the null model for the parametric bootstrap
##' @param parametric If TRUE a parametric bootstrap is calculated. If FALSE a
##' non-parametric (row-sampling) bootstrap is computed.
##' @param bollenstine Bollen-Stine transformation (non-parametric bootstrap) for bootstrap hypothesis testing.
##' @param constraints Logical indicating whether non-linear parameter
##' constraints should be included in the bootstrap procedure
##' @param sd Logical indicating whether standard error estimates should be
##' included in the bootstrap procedure
##' @param mc.cores Optional number of cores for parallel computing. If omitted future.apply will be used (see future::plan)
##' @param future.args arguments to future.apply::future_lapply
##' @param estimator String definining estimator, e.g. 'gaussian' (see
##' \code{estimator})
##' @param weights Optional weights matrix used by \code{estimator}
##' @param \dots Additional arguments, e.g. choice of estimator.
##' @aliases bootstrap.lvmfit
##' @usage
##'
##' \method{bootstrap}{lvm}(x,R=100,data,fun=NULL,control=list(),
##'                           p, parametric=FALSE, bollenstine=FALSE,
##'                           constraints=TRUE,sd=FALSE, mc.cores,
##'                           future.args=list(future.seed=TRUE),
##'                           ...)
##'
##' \method{bootstrap}{lvmfit}(x,R=100,data=model.frame(x),
##'                              control=list(start=coef(x)),
##'                              p=coef(x), parametric=FALSE, bollenstine=FALSE,
##'                              estimator=x$estimator,weights=Weights(x),...)
##'
##' @return A \code{bootstrap.lvm} object.
##' @author Klaus K. Holst
##' @seealso \code{\link{confint.lvmfit}}
##' @keywords models regression
##' @examples
##' m <- lvm(y~x)
##' d <- sim(m,100)
##' e <- estimate(lvm(y~x), data=d)
##' \donttest{ ## Reduce Ex.Timings
##' B <- bootstrap(e,R=50,mc.cores=1)
##' B
##' }
##' @export
bootstrap.lvm <- function(x, R = 100, data, fun = NULL, control = list(),
                          p, parametric = FALSE, bollenstine = FALSE,
                          constraints = TRUE, sd = FALSE, mc.cores,
                          future.args=list(future.seed=TRUE),
                          ...) {
  coefs <- sds <- c()
  on.exit(list(coef = coefs[-1, ], sd = sds[-1, ], coef0 = coefs[1, ], sd0 = sds[1, ], model = x))
  pb <- progressr::progressor(steps = R)
  pmis <- missing(p)
  bootfun <- function(i) {
    if (i == 0) {
      d0 <- data
    } else {
      if (!parametric | pmis) {
        d0 <- data[sample(seq_len(nrow(data)), replace = TRUE), ]
      } else {
        d0 <- sim(x, p = p, n = nrow(data))
      }
    }
    suppressWarnings(e0 <- estimate(x, data = d0, control = control, messages = 0, index = FALSE, ...))
    pb()
    if (!is.null(fun)) {
      coefs <- fun(e0)
      newsd <- NULL
    } else {
      coefs <- coef(e0)
      newsd <- c()
      if (sd) {
        newsd <- e0$coef[, 2]
      }
      if (constraints & length(constrain(x)) > 0) {
        cc <- constraints(e0, ...)
        coefs <- c(coefs, cc[, 1])
        names(coefs)[seq(length(coefs) - length(cc[, 1]) + 1, length(coefs))] <- rownames(cc)
        if (sd) {
          newsd <- c(newsd, cc[, 2])
        }
      }
    }
    return(list(coefs = coefs, sds = newsd))
  }
  if (bollenstine) {
    e0 <- estimate(x, data = data, control = control, messages = 0, index = FALSE, ...)
    mm <- modelVar(e0)
    mu <- mm$xi
    Y <- t(t(data[, manifest(e0)]) - as.vector(mu))
    Sigma <- mm$C
    S <- (ncol(Y) - 1) / ncol(Y) * var(Y)
    sSigma <- with(eigen(Sigma), vectors %*% diag(sqrt(values), ncol = ncol(vectors)) %*% t(vectors))
    isS <- with(eigen(S), vectors %*% diag(1 / sqrt(values), ncol = ncol(vectors)) %*% t(vectors))
    data <- as.matrix(Y) %*% (isS %*% sSigma)
    colnames(data) <- manifest(e0)
  }

  i <- 0
  if (!missing(mc.cores)) {
    res <- parallel::mclapply(0:R, bootfun, mc.cores=mc.cores)
  } else {
    res <- do.call(future_lapply, c(list(0:R, bootfun), future.args))
  }

  coefs <- matrix(unlist(lapply(res, function(x) x$coefs)), nrow = R + 1, byrow = TRUE)
  nn <- names(res[[1]]$coefs)
  if (!is.null(nn)) colnames(coefs) <- nn
  sds <- NULL
  if (sd) {
    sds <- matrix(unlist(lapply(res, function(x) x$sds)), nrow = R + 1, byrow = TRUE)
  }

  if (!is.null(fun)) {
    rownames(coefs) <- c()
    res <- list(coef = coefs[-1, , drop = FALSE], coef0 = coefs[1, ], model = x)
  } else {
    colnames(coefs) <- names(res[[1]]$coefs)
    rownames(coefs) <- c()
    if (sd) colnames(sds) <- colnames(coefs)
    res <- list(coef = coefs[-1, , drop = FALSE], sd = sds[-1, , drop = FALSE], coef0 = coefs[1, ], sd0 = sds[1, ], model = x, bollenstine = bollenstine)
  }
  class(res) <- "bootstrap.lvm"
  return(res)
}

##' @export
bootstrap.lvmfit <- function(x, R = 100, data = model.frame(x),
                             control = list(start = coef(x)),
                             p = coef(x), parametric = FALSE, bollenstine = FALSE,
                             estimator = x$estimator, weights = Weights(x), ...) {
  bootstrap.lvm(Model(x), R = R, data = data, control = control, estimator = estimator, weights = weights, parametric = parametric, bollenstine = bollenstine, p = p, ...)
}

##' @export
"print.bootstrap.lvm" <- function(x, idx, level = 0.95, ...) {
  cat("Non-parametric bootstrap statistics (R=", nrow(x$coef), "):\n\n", sep = "")
  uplow <- (c(0, 1) + c(1, -1) * (1 - level) / 2)
  nn <- paste(uplow * 100, "%")
  c1 <- t(apply(x$coef, 2, function(x) c(mean(x), sd(x), quantile(x, uplow))))

  c1 <- cbind(x$coef0, c1[, 1] - x$coef0, c1[, -1, drop = FALSE])
  colnames(c1) <- c("Estimate", "Bias", "Std.Err", nn)
  if (missing(idx)) {
    print(format(c1, ...), quote = FALSE)
  } else {
    print(format(c1[idx, , drop = FALSE], ...), quote = FALSE)
  }
  if (length(x$sd) > 0) {
    c2 <- t(apply(x$sd, 2, function(x) c(mean(x), sd(x), quantile(x, c(0.025, 0.975)))))
    c2 <- cbind(c2[, 1], c2[, 1] - x$sd0, c2[, -1])
    colnames(c2) <- c("Estimate", "Bias", "Std.Err", "2.5%", "97.5%")
    cat("\nStandard errors:\n")
    if (missing(idx)) {
      print(format(c2, ...), quote = FALSE)
    } else {
      print(format(c2[idx, , drop = FALSE], ...), quote = FALSE)
    }
  }
  cat("\n")
  invisible(x)
}
