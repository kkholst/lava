library(testthat)
library(lava)

set.seed(42)
n    <- 100
ic1  <- matrix(rnorm(n * 2), nrow=n, ncol=2)
ic2  <- matrix(rnorm(n * 2), nrow=n, ncol=2)
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
