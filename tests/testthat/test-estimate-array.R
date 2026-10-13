# estimate.default and estimate_calculus issues

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
  expect_equal(coef(e), expected, tolerance = 1e-8, ignore_attr = TRUE)
  x0 <- apply(x, 2, function(y) y-mean(y))
  var_expected <- crossprod(x0) / n^2
  expect_equal(vcov(e), var_expected, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("estimate.array variance formula is correct", {
  ## The variance computation in estimate.array currently evaluates
  ## mean(y - mean(y)^2) instead of mean((y - mean(y))^2).
  set.seed(1)
  x <- matrix(rnorm(200, mean = 3, sd = 2), ncol = 2)
  e <- estimate(x, type = "variance")
  expected <- apply(x, 2, function(y) mean((y - mean(y))^2))
  expect_equal(coef(e), expected, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("estimate.array variance handles NA without poisoning", {
  ## Inner mean(y) inside the variance formula does not pass na.rm=TRUE.
  set.seed(2)
  x <- matrix(rnorm(100), ncol = 2)
  x[1, 1] <- NA
  expect_no_error(e <- estimate(x, type = "variance"))
  expect_true(all(is.finite(coef(e))))
})

test_that("estimate.array quantile", {
  set.seed(1)
  y <- cbind(rnorm(100))
  q <- 0.75
  est <- quantile(y, type = 1L, probs=q)
  dens <- density(y)
  dens.est <- approxfun(dens)(est)
  infl <- (q - (y <= est)) / dens.est
  e1 <- estimate(coef = est, IC = infl)
  e2 <- estimate(y, type="quantile1", probs=q)
  # check manual vs quantile1 calc.
  expect_equal(IC(e1), IC(e2), ignore_attr = TRUE)
  expect_equal(coef(e1), coef(e2), ignore_attr = TRUE)
  # default for stats::quantile is type=7
  e3 <- estimate(y, type="quantile7", probs=q)
  expect_true(coef(e1) != coef(e3))
  expect_equal(coef(e3), quantile(y, probs=q), ignore_attr = TRUE)
  # check that multiple quantiles works
  qq <- c(0.25, 0.5, 0.75)
  e4 <- estimate(y, type="quantile7", probs=qq)
  expect_equal(coef(e4), quantile(y, qq), ignore_attr = TRUE)
})

test_that("estimate.array id", {
  set.seed(1)
  y <- cbind(rnorm(10))
  id <- 1:10
  e <- estimate(y, id=id, type="mean")
  expect_equal(index(e), id, ignore_attr = TRUE)
  # check correct order
  expect_equal(IC(e), y-mean(y), ignore_attr = TRUE)
  # string id
  ids <- rev(paste0("i", id))
  e <- estimate(y, id=ids)
  expect_equal(index(e), ids, ignore_attr = TRUE)
  expect_equal(IC(e), y-mean(y), ignore_attr = TRUE)
  # rownames id
  ys <- y
  rownames(ys) <- ids
  e <- estimate(ys)
  expect_equal(index(e), ids, ignore_attr = TRUE)
  expect_equal(IC(e), ys-mean(ys), ignore_attr = TRUE)
})
