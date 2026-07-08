context("merge.estimate")

library(testthat)
library(lava)

set.seed(42)
n    <- 100
ic1  <- matrix(rnorm(n * 2), nrow=n, ncol=2) |> scale(center=TRUE, scale=FALSE)
ic2  <- matrix(rnorm(n * 2), nrow=n, ncol=2) |> scale(center=TRUE, scale=FALSE)
e_ic1  <- estimate(coef=c(a=1, b=2), IC=ic1, id=1:n)
e_ic2  <- estimate(coef=c(c=3, d=4), IC=ic2, id=1:n+50)
e_noic  <- estimate(coef=c(e=5, f=6), vcov=diag(2))
e_noic2 <- estimate(coef=c(g=7, h=8), vcov=diag(2) * 2)

test_that("coefficients are concatenated correctly", {
  m <- merge(e_ic1, e_noic)
  expect_equal(coef(m), c(a=1, b=2, e=5, f=6))
})

test_that("vcov has correct dimensions", {
  m <- merge(e_ic1, e_noic)
  expect_equal(dim(vcov(m)), c(4L, 4L))
})

test_that("IC block is non-NA, cross-terms are NA", {
  m <- merge(e_ic1, e_noic)
  V <- vcov(m)
  expect_true(all(!is.na(V[1:2, 1:2])))
  expect_true(all(!is.na(V[3:4, 3:4])))
  expect_true(all(is.na(V[1:2, 3:4])))
  expect_true(all(is.na(V[3:4, 1:2])))
})

test_that("IC block matches direct merge of IC objects", {
  m_ic   <- merge(e_ic1, e_ic2)
  m_both <- merge(e_ic1, e_ic2, e_noic)
  expect_equal(vcov(m_both)[1:4, 1:4], vcov(m_ic))
})

test_that("non-IC block matches individual vcov", {
  m_both <- merge(e_ic1, e_ic2, e_noic)
  expect_equal(vcov(m_both)[5:6, 5:6], vcov(e_noic))
})

test_that("two non-IC objects: correct block placement and NA cross-terms", {
  m <- merge(e_ic1, e_noic, e_noic2)
  V <- vcov(m)
  expect_equal(dim(V), c(6L, 6L))
  expect_true(all(!is.na(V[1:2, 1:2])))
  expect_true(all(!is.na(V[3:4, 3:4])))
  expect_true(all(!is.na(V[5:6, 5:6])))
  expect_true(all(is.na(V[1:2, 3:6])))
  expect_true(all(is.na(V[3:4, c(1:2, 5:6)])))
  expect_true(all(is.na(V[5:6, 1:4])))
})

test_that("all non-IC objects: diagonal blocks filled, cross-terms NA", {
  m <- merge(e_noic, e_noic2)
  V <- vcov(m)
  expect_true(all(!is.na(V[1:2, 1:2])))
  expect_true(all(!is.na(V[3:4, 3:4])))
  expect_true(all(is.na(V[1:2, 3:4])))
  expect_true(all(is.na(V[3:4, 1:2])))
})

test_that("c.estimate with 'extra' attr", {
  e1 <- c(estimate(coef = c(a = 1), vcov = matrix(0.1)), e1.iter = 2)
  e2 <- c(estimate(coef = c(b = 2), vcov = matrix(0.2)), e2.iter = 3)
  e3 <- c(estimate(coef = c(c = 2), vcov = matrix(0.3)))

  e <- c(e1, e2) # two estimate.extra objects
  expect_equal(c(coef(e1), coef(e2)), coef(e))
  expect_equal(c(vcov(e1), vcov(e2)),
               diag(vcov(e)), check.attributes=FALSE)
  expect_equal(c(attr(e1, "extra"), attr(e2, "extra")),
                 attr(e, "extra"))

  e <- c(e1, e3) # estimatee.xtra + estimate
  expect_equal(c(coef(e1), coef(e3)), coef(e))
  expect_equal(c(vcov(e1), vcov(e3)),
               diag(vcov(e)), check.attributes=FALSE)
  expect_equal(attr(e1, "extra"), attr(e, "extra"))

  e <- c(e3, e1) # estimate + estimate.extra
  expect_equal(c(coef(e3), coef(e1)), coef(e))
  expect_equal(c(vcov(e3), vcov(e1)),
               diag(vcov(e)), check.attributes=FALSE)
  expect_equal(attr(e1, "extra"), attr(e, "extra"))

  e <- c(e1, e2, d=2) # estimate + scalar/extra
  expect_equal(c(coef(e1), coef(e2)), coef(e))
  expect_equal(c(vcov(e1), vcov(e2)),
               diag(vcov(e)), check.attributes=FALSE)
  expect_equal(c(attr(e1, "extra"), attr(e2, "extra"), d=2),
               attr(e, "extra"))

  e <- c(e3, c(a = 1, b = 2)) # named scalar vector via unnamed argument
  expect_equal(c(a = 1, b = 2), attr(e, "extra"))

  e <- c(e3, e1.iter = 1) # duplicated name of extra variable
  e <- c(e, e1)
  expect_equal(c(e1.iter = 1, e1.iter = 2), attr(e, "extra"))

  e <- c(e1, converged = FALSE) # boolean extra variable
  # internal conversion to numeric
  expect_equal(c(e1.iter = 2, converged = 0), attr(e, "extra"))
})

test_that("estimate preserves original order with stack=TRUE", {
  set.seed(1)
  n <- 10
  x <- seq_len(n)
  y <- rnorm(n, sd=0.1) + x
  id1 <- as.character(x + 5)
  d <- data.frame(y, x)
  rownames(d) <- id1

  g <- glm(y~x, data=d)
  e <- estimate(g, id=id1)
  X <- model.matrix(g)
  r <- residuals(g)
  ic1 <- cbind(r, r*d$x) %*% solve(crossprod(X))*nrow(X)
  ic2 <- IC(e)
  expect_equivalent(ic1,ic2)
})

test_that("merge preserves order (of first IC)", {
  set.seed(1)
  id1 <- as.character(1:10 + 5)
  id2 <- as.character(1:10)
  a <- estimate(coef=1, IC=scale(rnorm(length(id1))), id=id1)
  b <- estimate(coef=2, IC=scale(rnorm(length(id2))), id=id2)

  expect_equivalent(id1, rownames(IC(a)))
  expect_equivalent(id2, rownames(IC(b)))

  e <- merge(a, b)
  expect_equivalent(rownames(IC(e))[seq_along(id1)],
                    rownames(IC(a)))

  ic1 <- IC(e)[seq_along(id1),1,drop=FALSE]
  expect_equivalent(IC(a)* 3/2, ic1) # same expect for IPW due to incomplete data
})
