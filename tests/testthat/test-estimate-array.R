context("estimate.default and estimate_calculus issues")

library(testthat)
library(lava)


test_that("estimate.array mean formula is correct", {
  ## The variance computation in estimate.array currently evaluates
  ## mean(y - mean(y)^2) instead of mean((y - mean(y))^2).
  set.seed(1)
  n <- 100
  x <- matrix(rnorm(2*n, mean = 3, sd = 2), ncol = 2)
  e <- estimate(x)
  expected <- apply(x, 2, function(y) mean(y))
  expect_equivalent(coef(e), expected, tolerance = 1e-8)
  x0 <- apply(x, 2, function(y) y-mean(y))
  var_expected <- crossprod(x0) / n^2
  expect_equivalent(vcov(e), var_expected, tolerance = 1e-8)
})

test_that("estimate.array variance formula is correct", {
  ## The variance computation in estimate.array currently evaluates
  ## mean(y - mean(y)^2) instead of mean((y - mean(y))^2).
  set.seed(1)
  x <- matrix(rnorm(200, mean = 3, sd = 2), ncol = 2)
  e <- estimate(x, type = "variance")
  expected <- apply(x, 2, function(y) mean((y - mean(y))^2))
  expect_equivalent(coef(e), expected, tolerance = 1e-8)
})

test_that("estimate.array variance handles NA without poisoning", {
  ## Inner mean(y) inside the variance formula does not pass na.rm=TRUE.
  set.seed(2)
  x <- matrix(rnorm(100), ncol = 2)
  x[1, 1] <- NA
  expect_no_error(e <- estimate(x, type = "variance"))
  expect_true(all(is.finite(coef(e))))
})
