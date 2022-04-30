##' Model searching
##'
##' Performs Wald or score tests
##'
##'
##' @aliases modelsearch
##' @param x \code{lvmfit}-object
##' @param k Number of parameters to test simultaneously. For \code{equivalence}
##' the number of additional associations to be added instead of \code{rel}.
##' @param dir Direction to do model search. "forward" := add
##' associations/arrows to model/graph (score tests), "backward" := remove
##' associations/arrows from model/graph (wald test)
##' @param type If equal to 'correlation' only consider score tests for covariance parameters. If equal to 'regression' go through direct effects only  (default 'all' is to do both)
##' @param ... Additional arguments to be passed to the low level functions
##' @return Matrix of test-statistics and p-values
##' @author Klaus K. Holst
##' @seealso \code{\link{compare}}, \code{\link{equivalence}}
##' @keywords htest
##' @examples
##'
##' m <- lvm();
##' regression(m) <- c(y1,y2,y3) ~ eta; latent(m) <- ~eta
##' regression(m) <- eta ~ x
##' m0 <- m; regression(m0) <- y2 ~ x
##' dd <- sim(m0,100)[,manifest(m0)]
##' e <- estimate(m,dd);
##' modelsearch(e,messages=0)
##' modelsearch(e,messages=0,type="cor")
##' @export
modelsearch <- function(x, k = 1, dir = "forward", type = "all", ...) {
  if (dir == "forward") {
    res <- forwardsearch(x, k, type = type, ...)
    return(res)
  }
  if (dir == "backstep") {
    res <- backwardeliminate(x, ...)
    return(res)
  }
  res <- backwardsearch(x, k, ...)
  return(res)
}

backwardeliminate <- function(x,
                              keep = NULL,
                              pthres = 0.05,
                              AIC = FALSE,
                              messages = 0,
                              missing = FALSE,
                              intercepts = FALSE,
                              maxsteps = Inf,
                              information = "E",
                              data,
                              ...) {
  if (inherits(x, "lvm")) {
    M <- x
  } else {
    M <- Model(x)
  }
  if (missing(data)) data <- model.frame(x)

  dots <- list(...)
  if (is.null(dots$control$start)) {
    p0 <- estimate(M, data, quick = TRUE, messages = messages, missing = FALSE, ...)
    dots$control <- c(dots$control, list(start = p0, information = "E"))
  }

  ff <- function() {
    ii <- grep("m", names(coef(M)))
    vv <- variances(M, mean = TRUE)
    args <- c(list(x = M, data = data, missing = missing, quick = TRUE, messages = messages), dots)
    cc <- do.call("estimate", args)
    if (is.numeric(cc)) {
      I0 <- information(M, p = cc, data = data, type = information)[-c(ii, vv), -c(ii, vv)]
      cc0 <- cc[-c(ii, vv)]
      res <- (pnorm(abs(cc0 / sqrt(diag(solve(I0)))), lower.tail = FALSE)) * 2
      attributes(res)$coef <- cc
    } else {
      coefs <- coef(cc)
      res <- (pnorm(abs(coefs / sqrt(diag(vcov(cc)))), lower.tail = FALSE)) * 2
      res <- res[-c(ii, vv)]
      attributes(res)$coef <- coefs
    }
    return(res)
  }

  done <- FALSE
  i <- 0
  while (!done & i < maxsteps) {
    p <- ff()
    ordp <- order(p, decreasing = TRUE)
    curp <- p[ordp[1]]
    if (curp < pthres) break
    dots$control$start <- attributes(p)$coef[-ordp[1]]
    if (messages) message("Removed: ", names(curp), " p-value: ", round(curp, 3))
    ## var1 <- unlist(strsplit(names(curp),lava.options()$symbol[1]))
    nn <- strsplit(names(curp), paste0(lava.options()$symbol, collapse = "|"))[[1]]
    cancel(M) <- nn
  }

  if (messages) message("")
  return(M)
}

backwardsearch <- function(x, k = 1, ...) {
  if (!inherits(x, "lvmfit")) stop("Expected an object of class 'lvmfit'.")
  p <- pars(x)
  cur <- Model(x)
  pp <- modelPar(cur, p)
  p1 <- pp$p
  Tests <- c()
  Vars <- list()

  parnotvar <- setdiff(seq_along(p1), variances(Model(x))) ## We don't want to perform tests on the boundary of the parameter space
  freecomb <- utils::combn(parnotvar, k)

  for (i in seq_len(ncol(freecomb)))
  {
    cc0 <- coef(cur, mean = FALSE, messages = 0, symbol = lava.options()$symbol)
    ii <- freecomb[, i]
    p0 <- p1
    p0[ii] <- 0
    R <- diag(nrow = length(p0))
    R <- matrix(R[ii, ], nrow = length(ii))
    I <- information(Model(x), p = p1, n = x$data$n, data = model.frame(x))
    if (!is.null(pp$meanpar)) {
      rmidx <- seq_along(pp$meanpar)
      I <- I[-rmidx, -rmidx]
    }
    iI <- solve(I)
    W <- t(rbind(R) %*% p1) %*% solve(R %*% iI %*% t(R)) %*% (cbind(R) %*% p1)
    Tests <- c(Tests, W)
    Vars <- c(Vars, list(cc0[ii]))
  }
  ord <- order(Tests, decreasing = TRUE)
  Tests <- cbind(Tests, pchisq(Tests, k, lower.tail = FALSE))
  colnames(Tests) <- c("Test Statistic", "P-value")
  res <- list(test = Tests[ord, , drop = FALSE], var = Vars[ord])
  PM <- matrix(ncol = 3, nrow = 0)
  for (i in seq_len(nrow(Tests))) {
    if (!is.na(res$test[i, 1])) {
      newrow <- c(formatC(res$test[i, 1]), formatC(res$test[i, 2]), paste(res$var[[i]], collapse = ", "))
      PM <- rbind(PM, newrow)
    }
  }
  colnames(PM) <- c("Wald: W", "P(W>w)", "Index")
  rownames(PM) <- rep("", nrow(PM))

  res <- list(res = PM, test = res$test)
  class(res) <- "modelsearch"
  res
}

forwardsearch <- function(x, k = 1, messages = lava.options()$messages, type = "all", exclude.var = NULL, ...) {
  if (!inherits(x, "lvmfit")) stop("Expected an object of class 'lvmfit'.")

  p <- pars(x, reorder = TRUE)
  cur <- Model(x)
  Y <- endogenous(x)
  X <- exogenous(x)
  V <- vars(x)
  q <- length(Y)
  qx <- length(X)
  npar.sat <- q + q * (q - 1) / 2 + q * qx
  npar.cur <- index(cur)$npar
  nfree <- npar.sat - npar.cur
  if (nfree < k) {
    message("Cannot free ", k, " variables from model.\n")
    return()
  }

  directional <- !(tolower(type) %in% c("cor", "correlation", "cov", "covariance"))
  all <- tolower(type) %in% c("all", "both")

  Tests <- c()
  Vars <- list()
  AP <- with(index(cur), A + t(A) + P)
  restricted <- c()
  ## idx1 <- seq_len(ncol(AP)-1)
  ## idx2 <- seq(i+1,nrow(AP))

  idx <- seq_len(ncol(AP))
  if (!is.null(exclude.var)) {
    if (is.character(exclude.var)) {
      exclude.var <- match(exclude.var, V)
    }
    idx <- setdiff(idx, exclude.var)
  }
  for (i0 in seq_len(length(idx) - 1)) {
    for (j0 in seq(i0 + 1, length(idx))) {
      i <- idx[i0]
      j <- idx[j0]
      if (AP[j, i] == 0) {
        restricted <- rbind(restricted, c(i, j))
      }
    }
  }

  if (is.null(restricted)) {
    return(NULL)
  }
  if (all) {
    ntest <- nrow(restricted)
    directional <- rep(FALSE, ntest)
    restricted <- rbind(restricted, restricted, restricted[, 2:1])
    directional <- c(directional, rep(TRUE, 2 * ntest))
  } else {
    if (directional) {
      restricted <- rbind(restricted, restricted[, 2:1])
    }
    directional <- rep(directional, nrow(restricted))
  }

  restrictedcomb <- utils::combn(seq_len(nrow(restricted)), k) # Combinations of k-additions to the model

  if (!inherits(model.frame(x), c("data.frame", "matrix"))) {
    n <- model.frame(x)$n
    S <- model.frame(x)$S
    mu <- model.frame(x)$mu
  } else {
    n <- nrow(model.frame(x))
    S <- (n - 1) / n * var(model.frame(x), na.rm = TRUE)
    mu <- colMeans(model.frame(x), na.rm = TRUE)
  }
  pb <- progressr::progressor(steps = ncol(restrictedcomb))
  if (messages > 0) {
    message("Calculating score test for ", ncol(restrictedcomb), " models:")
  }
  for (i in seq_len(ncol(restrictedcomb))) {
    if (messages > 0) pb()
    varlist <- c()
    altmodel <- cur ## HA: altmodel, H0: cur
    for (j in seq_len(k)) {
      myvar <- restricted[restrictedcomb[j, i], ]
      if (any(wx <- V[myvar] %in% X)) {
        altmodel <- regression(altmodel, V[myvar][which(!wx)], V[myvar][which(wx)])
      } else {
        if (directional[i]) {
          covariance(altmodel, pairwise = TRUE) <- V[myvar]
        }
        covariance(altmodel, pairwise = TRUE) <- V[myvar]
      }
      varlist <- rbind(varlist, V[myvar])
    }
    altmodel$parpos <- NULL
    altmodel <- updatelvm(altmodel, deriv = TRUE, zeroones = TRUE, mean = TRUE)
    cc <- coef(altmodel, mean = TRUE, messages = 0, symbol = lava.options()$symbol)
    cc0 <- coef(cur, mean = TRUE, messages = 0, symbol = lava.options()$symbol)
    p1 <- numeric(length(p) + k)
    ## Need to be sure we place 0 at the correct position
    for (ic in seq_along(cc)) {
      idx <- match(cc[ic], cc0)
      if (!is.na(idx)) {
        p1[ic] <- p[idx]
      }
    }
    if (x$estimator == "gaussian" && !inherits(x, "lvm.missing")) {
      Sc2 <- score(altmodel,
        p = p1, data = NULL,
        model = x$estimator, weights = Weights(x), S = S, mu = mu, n = n
      )
    } else {
      Sc2 <- score(altmodel,
        p = p1, data = model.frame(x),
        model = x$estimator, weights = Weights(x)
      )
    }
    I <- information(altmodel, p = p1, n = n, data = model.frame(x), weights = Weights(x), estimator = x$estimator) ## [-rmidx,-rmidx]

    iI <- try(Inverse(I), silent = TRUE)
    Q <- ifelse(inherits(iI, "try-error"), NA, ## Score test
      (Sc2) %*% iI %*% t(Sc2)
    )
    Tests <- c(Tests, Q)
    Vars <- c(Vars, list(varlist))
  }

  ##if (messages > 0) close(pb)
  ord <- order(Tests)
  Tests <- cbind(Tests, pchisq(Tests, k, lower.tail = FALSE))
  colnames(Tests) <- c("Test Statistic", "P-value")
  PM <- c()
  for (i in seq_len(nrow(Tests))) {
    if (!is.na(Tests[i, 1])) {
      vv <- apply(Vars[[i]], 1, function(x) paste(x, collapse = lava.options()$symbol[2 - directional[i]]))
      newrow <- c(
        formatC(Tests[i, 1]), formatC(Tests[i, 2]),
        paste(vv, collapse = ",")
      )
      PM <- rbind(PM, newrow)
    }
  }
  if (is.null(PM)) {
    message("Saturated model")
    return(invisible(NULL))
  }
  Tests <- Tests[ord, , drop = FALSE]
  Vars <- Vars[ord]
  PM <- PM[ord, , drop = FALSE]

  colnames(PM) <- c("Score: S", "P(S>s)", "Index")
  rownames(PM) <- rep("", nrow(PM))
  res <- list(res = PM, test = Tests, var = Vars, directional = directional)
  class(res) <- "modelsearch"
  return(res)
}

##' @export
print.modelsearch <- function(x, tail = nrow(x$res), adj = c("holm", "BH"), ...) {
  N <- nrow(x$res)
  if (!is.null(adj)) {
    ##    adjp <- rev(holm(as.numeric(x$test[,2])))
    adjp <- rbind(sapply(adj, function(i) p.adjust(x$test[, 2], method = i)))
    colnames(adjp) <- adj
    x$res <- cbind(x$res, rbind(formatC(adjp)))
  }
  print(x$res[seq(N - tail + 1, N), ], quote = FALSE, ...)
  invisible(x)
}
