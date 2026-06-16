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

test_That("c.estimate.extra", {
  e1 <- c(estimate(coef = c(a = 1), vcov = matrix(0.1)), e1.iter = 2)
  e2 <- c(estimate(coef = c(b = 2), vcov = matrix(0.2)), e2.iter = 3)
  e3 <- c(estimate(coef = c(c = 2), vcov = matrix(0.3)))

  e <- c(e1, e2) # two estimate.extra objects
  expect_equal(c(coef(e1$estimate), coef(e2$estimate)), coef(e$estimate))
  expect_equal(c(vcov(e1$estimate), vcov(e2$estimate)),
               diag(vcov(e$estimate)), check.attributes=FALSE)
  expect_equal(c(e1$extra, e2$extra), e$extra)

  e <- c(e1, e3) # estimate.extra + estimate
  expect_equal(c(coef(e1$estimate), coef(e3)), coef(e$estimate))
  expect_equal(c(vcov(e1$estimate), vcov(e3)),
               diag(vcov(e$estimate)), check.attributes=FALSE)
  expect_equal(c(e1$extra), e$extra)

  e <- c(e3, e1) # estimate + estimate.extra
  expect_equal(c(coef(e3), coef(e1$estimate)), coef(e$estimate))
  expect_equal(c(vcov(e3), vcov(e1$estimate)),
               diag(vcov(e$estimate)), check.attributes=FALSE)
  expect_equal(c(e1$extra), e$extra)

  e <- c(e1, e2, d=2) # estimate + scalar/extra
  expect_equal(c(coef(e1$estimate), coef(e2$estimate)), coef(e$estimate))
  expect_equal(c(vcov(e1$estimate), vcov(e2$estimate)),
               diag(vcov(e$estimate)), check.attributes=FALSE)
  expect_equal(c(e1$extra, e2$extra, d=2), e$extra)
})
