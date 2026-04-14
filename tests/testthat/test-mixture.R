library(testthat)

context("mixture lvm")

setup_mixture_data <- function(seed = 42, n = 500) {
  set.seed(seed)
  m0 <- lvm(list(y ~ x + z, x ~ z))
  distribution(m0, ~z) <- binomial.lvm()
  d <- sim(m0, n, p = c("y~z" = 2, "y~x" = 1))
  list(data = d, model = lvm(y ~ x, x ~ 1))
}

s <- setup_mixture_data()
m <- baptize(s$model)
intercept(m, ~x + y) <- NA
result <- mixture(m, k = 2, data = s$data,
                    control = list(tol = 1e-4, trace = 0))

test_that("mixture: result contains expected components", {
  expect_s3_class(result, "lvm.mixture")
  expect_true(!is.null(result$theta))
  expect_true(!is.null(result$prob))
  expect_true(!is.null(result$gamma))
  expect_true(!is.null(result$member))
  expect_true(!is.null(result$vcov))
  expect_true(!is.null(result$k))
  ## mixture: correct number of mixture components (k=2)
  expect_equal(result$k, 2)
})

test_that("mixture: mixing probabilities sum to 1", {
  prob <- tail(result$prob, 1)  # last row = final estimates
  expect_equal(sum(prob), 1, tolerance = 1e-6)
  expect_true(all(prob >= 0))
  expect_true(all(prob <= 1))
})

test_that("mixture: posterior probabilities (gamma)", {
  #  sum to 1 per observation
  row_sums <- rowSums(result$gamma)
  expect_equal(row_sums, rep(1, nrow(s$data)), tolerance = 1e-6)
  ## gamma has correct dimensions
  expect_equal(dim(result$gamma), c(nrow(s$data), 2))
})

test_that("coef.mixture", {
  cc <- coef(result)
  expect_true(is.numeric(cc))
  expect_true(!is.null(names(cc)))
  ## coef(list=TRUE) returns a list of length k
  cc_list <- coef(result, list = TRUE)
  expect_true(is.list(cc_list))
  expect_equal(length(cc_list), 2)
})

test_that("logLik.mixture", {
  ll <- logLik(result)
  expect_true(is.finite(ll))
  expect_equal(class(ll), "logLik")
  expect_equal(attr(ll, "nobs"), nrow(s$data))
})

test_that("vcov.mixture", {
  V <- vcov(result)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), ncol(V))
  np <- length(coef(result))
  expect_equal(dim(vcov(result)), c(np, np))
})

test_that("mixture: score() at MLE is approximately zero", {
  sc <- score(result)
  expect_true(sum(sc^2) < 1)  # score norm should be small at MLE
  sc_indiv <- score(result, indiv = TRUE)
  np <- length(coef(result))
  expect_equal(dim(sc_indiv), c(nrow(s$data), np))
})

test_that("mixture: summary() returns object of class 'summary.lvm.mixture'", {
  s_result <- summary(result)
  expect_s3_class(s_result, "summary.lvm.mixture")
})

test_that("mixture: summary() AIC is finite", {
  s_result <- summary(result)
  expect_true(is.finite(s_result$AIC))
})
