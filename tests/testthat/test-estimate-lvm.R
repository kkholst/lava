context("estimate lvm")

test_that("estimate.lvm cluster", {
  set.seed(1)
  m <- lvm(y ~ x)
  d <- sim(m,100)
  l <- lm(y ~ x, data=d)
  e <- estimate(m, data=d, id=1:nrow(d))
  e2 <- estimate(l, id=1:nrow(d))
  expect_true(mean((vcov(e2)-vcov(e)[c(1,2),c(1,2)])^2)<1e-6)
})

test_that("estimate.list", {
  d <- data.frame(y=rnorm(10), x=rnorm(10))
  l <- lm(y~x, data=d)
  e <- estimate(list(l, l, l))
  cc <- coef(e)
  expect_true(length(cc)==3L)
  expect_equivalent(cc[[1]], coef(l))
  expect_true(inherits(e[[1]], "estimate"))
})
