pzmax <- function(alpha, S) {
  ## P(Zmax > z) Family wise error rate, Zmax = max |Z_i|
  if (!requireNamespace("mets",quietly=TRUE))
    stop("'mets' package required")
  k <- nrow(S)
  zz <- qnorm(1-alpha/2)
  unlist(lapply(zz, function(z)
    1 - mets::pmvn(lower=rep(-z, k),
                   upper=rep(z, k),
                   sigma=cov2cor(S))
    ))
}

##' @export
alpha_zmax <- function(object, method, alpha = 0.05, ...) {
  if (!inherits(object, "estimate")) {
    stop("Expected an 'estimate' object")
  }
  dots <- list(...)
  null <- object$compare$null
  if (!("null" %in% names(dots)) && !is.null(null)) {
    dots$null <- null
  }
  object <- do.call(estimate, c(list(object), dots))
  est <- parameter(object)[, c(1, 5), drop = FALSE]
  padj <- pzmax(est[, 2], vcov(object))
  res <- cbind(est, padj)
  colnames(res)[3] <- paste0("Adj.", colnames(res)[2])
  f <- function(a) pzmax(a, vcov(object)) - alpha
  adj <- uniroot(f, lower = 0, upper = alpha)$root
  attributes(res)["adjusted.significance.level"] <- adj
  res
}

##' Closed testing procedure
##'
##' Given p hypotheses H1, ..., Hp all 2^p-1 intersection hypotheses are
##' calculated and adjusted p-values are obtained for Hj is calculated as the
##' max p-value of all intersection hypotheses containing Hj. Example, for p=3,
##' the adjusted p-value for H1 will be obtained from \{(H1, H2, H3), (H1,H2),
##' (H1,H3), (H1)\}.
##' @aliases closed_testing alpha_zmax
##' @param object `estimate` object
##' @param test function that conducts hypothesis test. See details below.
##' @param ... Additional arguments passed to `test`
##' @export
##' @details The function `test` should be a function `function(object, index,
##'   ...)` which as its first argument takes an `estimate` object and and wit
##'   an argument `index` which is a integer vector specifying which
##'   subcomponents of `object` to test. The ellipsis argument can be any other
##'   arguments used in the test function. The function \code{test_wald} is an
##'   example of valid test function (which has an additional argument `null` in
##'   reference to the above mentioned ellipsis arguments).
##' @references Marcus, R; Peritz, E; Gabriel, KR (1976).
##'   "On closed testing procedures with special reference to ordered analysis
##'   of variance". Biometrika. 63 (3): 655–660.
##' @examples
##' m <- lvm()
##' regression(m, c(y1,y2,y3,y4)~x) <- c(0, 0.25, 0, 0.25)
##' regression(m, to=endogenous(m), from="u") <- 1
##' variance(m,endogenous(m)) <- 1
##' set.seed(1)
##' d <- sim(m, 200)
##' l1 <- lm(y1~x,d)
##' l2 <- lm(y2~x,d)
##' l3 <- lm(y3~x,d)
##' l4 <- lm(y4~x,d)
##'
##' (a <- merge(l1, l2, l3, l4, subset=2))
##' if (requireNamespace("mets",quietly=TRUE)) {
##'    alpha_zmax(a)
##' }
##' adj <- closed_testing(a)
##' adj
##' adj$p.value
##' summary(adj)
closed_testing <- function(object, test = test_wald, ...) {
  if (!inherits(object, "estimate")) stop("`estimate` object needed")
  idx <- seq_along(coef(object))
  if (length(idx) > 15)
    stop("Too many tests. Consider some other adjustment method")
  combs <- pvals <- c()
  for (i in seq_along(idx)) {
    co <- combn(idx, i)
    pp <- numeric(ncol(co))
    for (j in seq_along(pp)) {
      ii <- co[, j]
      tt <- test(object, index = ii, ...)
      pp[j] <- tt$p.value
    }
    combs <- c(combs, list(co))
    pvals <- c(pvals, list(pp))
  }
  pmax <- c()
  for (k in seq_along(idx)) {
    pk <- c()
    for (i in seq_along(idx)) {
      cols <- apply(combs[[i]], 2, function(x) k %in% x)
      pk <- c(pk, pvals[[i]][which(cols)])
    }
    pmax <- c(pmax, max(pk))
  }
  res <- cbind(coef(object), pmax)
  colnames(res) <- c("Estimate", "adj.p")
  rownames(res) <- names(coef(object))
  names(pmax) <- names(coef(object))
  structure(
    list(
      adjp = res,
      call = match.call(),
      args = list(...),
      hypotheses = combs,
      raw.pval = pvals,
      estimate = coef(object),
      p.value = pmax),
    class = "test_adj"
  )
}


#' @export
print.test_adj <- function(x, ...) {
  cat("Call: ")
  print(x$call, quote = FALSE)
  cat("\n")
  print(x$adjp)
}

#' @export
summary.test_adj <- function(object, ...) {
  cat("Call: ")
  print(object$call, quote = FALSE)
  # Print adjusted p-values
  cli::cli_h2("Adjusted p-values")
  print(object$adjp)
  cat("\n")
  # Print raw p-values for intersection hypotheses
  cli::cli_h2("Raw p-values for intersection hypotheses")
  param_names <- names(object$estimate)
  for (i in seq_along(object$hypotheses)) {
    cat(sprintf("%d-way intersections:\n", i))
    comb_matrix <- object$hypotheses[[i]]
    pvals <- object$raw.pval[[i]]
    for (j in seq_len(ncol(comb_matrix))) {
      indices <- comb_matrix[, j]
      hyp_names <- param_names[indices]
      # Create hypothesis notation
      hyp_str <- paste0("{", paste(hyp_names, collapse = ", "), "}")
      # Print with p-value
      cat(sprintf("  %-40s p = %.4f\n", hyp_str, pvals[j]))
    }
    cat("\n")
  }
  invisible(object)
}

#' @export
test_wald <- function(par,
                      vcov,
                      null = NULL,
                      index = NULL) {
  if (!inherits(par, "estimate")) {
    par <- lava::estimate(par=par, vcov=vcov)
  }
  if (is.null(null)) {
    null <- rep(0, length(coef(par)))
  }
  B <- diag(length(coef(par)))
  if (!is.null(index)) {
    if (length(coef(par)) < length(index)) stop("wrong `index`")
    null <- null[index]
    B <- B[index, , drop=FALSE]
  }
  lava::compare(par, contrast = B, null=null)
}
