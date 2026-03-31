context("Inference")

set.seed(1)
a1 <- estimate(coef=1, IC=rnorm(10), id=1:10, labels="a1")
a2 <- estimate(coef=2, IC=rnorm(10), id=1:10, labels="a2")
a <- merge(a1, a2)

b1 <- estimate(coef=0.5, IC=rnorm(10), id=1:10, labels="b1")
b2 <- estimate(coef=0.9, IC=rnorm(10), id=1:10, labels="b2")
b <- merge(b1, b2)

test_that("+.estimate", {
  ## check addition: estimate + numeric with dim(estimate)=1
  e1 <- a1 + c(2,4,6)
  e2 <- estimate(a1, \(p) p + c(2,4,6))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  ## check addition: estimate + numeric with dim(estimate)>1
  e1 <- a + 2
  e2 <- estimate(a, \(p) p + 2)
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- a + c(2,4)
  e2 <- estimate(a, \(p) p + c(2,4))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  # commutative
  e1 <- c(2,4) + a # 2.dim. estimate
  e2 <- estimate(a, \(p) p + c(2,4))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- c(2,4) + a1 # 1.dim
  e2 <- estimate(a1, \(p) p + c(2,4))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  # addition with two estimate objects
  e1 <- a1 + a2
  e2 <- estimate(a, cbind(1, 1))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- a + b # 2d. estimates
  e2 <- estimate(merge(a, b), rbind(c(1, 0, 1, 0), c(0, 1, 0, 1)))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
})

test_that("-.estimate", {
  ## check minus operator
  e1 <- -a
  e2 <- estimate(a, function(p) -p)
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  ## estimate + numeric with dim(estimate)=1
  e1 <- a1 - c(2,4,6)
  e2 <- estimate(a1, \(p) p - c(2,4,6))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  ## estimate + numeric with dim(estimate)>1
  e1 <- a - 2
  e2 <- estimate(a, \(p) p - 2)
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- a - c(2,4)
  e2 <- estimate(a, \(p) p - c(2,4))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  # reverse order of numeric and estimate
  e1 <- c(2,4) - a # 2.dim. estimate
  e2 <- estimate(a, \(p) c(2,4) - p)
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- c(2,4) - a1 # 1.dim
  e2 <- estimate(a1, \(p) c(2,4) - p)
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  # subtraction of two estimate objects
  e1 <- a1 - a2
  e2 <- estimate(a, cbind(1, -1))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- a - b # 2d. estimates
  e2 <- estimate(merge(a, b), rbind(c(1, 0, -1, 0), c(0, 1, 0, -1)))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
})

test_that("*.estimate", {
  ## check multiplication: estimate + numeric with dim(estimate)=1
  e1 <- a1 * c(2,4,6)
  e2 <- estimate(a1, \(p) p * c(2,4,6))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  ## estimate + numeric with dim(estimate)>1
  e1 <- a * 2
  e2 <- estimate(a, \(p) p * 2)
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- a * c(2,4)
  e2 <- estimate(a, \(p) p * c(2,4))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  # change order
  e1 <- c(2,4) * a # 2.dim. estimate
  e2 <- estimate(a, \(p) p * c(2,4))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- c(2,4) * a1 # 1.dim
  e2 <- estimate(a1, \(p) p * c(2,4))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  # with two estimate objects
  e1 <- a1 * a2
  e2 <- estimate(a, prod)
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- a * b # 2d. estimates
  e2 <- estimate(merge(a, b), function(p) c(p[1]*p[3], p[2]*p[4]))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
})

test_that("/.estimate", {
  ## check division: estimate + numeric with dim(estimate)=1
  e1 <- a1 / c(2,4,6)
  e2 <- estimate(a1, \(p) p / c(2,4,6))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  ## estimate + numeric with dim(estimate)>1
  e1 <- a / 2
  e2 <- estimate(a, \(p) p / 2)
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- a / c(2,4)
  e2 <- estimate(a, \(p) p / c(2,4))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  # numeric first
  e1 <- c(2,4) / a # 2.dim. estimate
  e2 <- estimate(a, \(p) c(2,4) / p)
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- c(2,4) / a1 # 1.dim
  e2 <- estimate(a1, \(p) c(2,4) / p)
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  # with two estimate objects
  e1 <- a1 / a2
  e2 <- estimate(a, function(p) p[1] / p[2])
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- a / b # 2d. estimates
  e2 <- estimate(merge(a, b), function(p) c(p[1]/p[3], p[2]/p[4]))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
})

test_that("^.estimate", {
  ## check power-func: estimate + numeric with dim(estimate)=1
  e1 <- a1^0.25
  e2 <- estimate(a1, function(p) p**0.25)
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- a1 / c(2,4,6)
  e2 <- estimate(a1, \(p) p / c(2,4,6))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  ## estimate + numeric with dim(estimate)>1
  e1 <- a ^ 3
  e2 <- estimate(a, \(p) p ^ 3)
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- a ^ c(2,4)
  e2 <- estimate(a, \(p) p ^ c(2,4))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  # numeric first
  e1 <- c(2,4) ^ a # 2.dim. estimate
  e2 <- estimate(a, \(p) c(2,4) ^ p)
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- c(2,4) ^ a1 # 1.dim
  e2 <- estimate(a1, \(p) c(2,4) ^ p)
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  # addition with two estimate objects
  e1 <- a1 ^ a2
  e2 <- estimate(a, function(p) p[1] ^ p[2])
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- a ^ b # 2d. estimates
  e2 <- estimate(merge(a, b), function(p) c(p[1]^p[3], p[2]^p[4]))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
})

test_that("math functions", {
  e <- c(a, b) # merge
  e2 <- merge(a, b); e2$model.index <- NULL
  expect_equivalent(e, e2)

  # subset and merge
  e1 <- c(a=a["a1"], b=b["b1"])
  expect_true(length(coef(e1)) == 2L)
  expect_equivalent(names(coef(e1)), c("a", "b"))

  # product
  e1 <- prod(e, labels="prod")
  e2 <- estimate(e, function(p) prod(p))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  # sum
  e1 <- sum(e, labels="sum")
  e2 <- estimate(e, function(p) sum(p))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  # inner-prod
  e1 <- a %*% b
  e2 <- estimate(e, function(p) sum(p[1:2]*p[3:4]))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- a %*% c(2, 3)
  e2 <- estimate(a, function(p) sum(p[1]*2 + p[2]*3))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  # log
  e1 <- log(e["a1"])
  e2 <- estimate(e, function(p) log(p[1]))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- log(e)
  e2 <- estimate(e, function(p) log(p))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  e1 <- log1p(b)
  e2 <- estimate(b, function(p) log1p(p))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  # exp
  e1 <- exp(e["a1"])
  e2 <- estimate(e, function(p) exp(p[1]))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
  e1 <- exp(e)
  e2 <- estimate(e, function(p) exp(p))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  e1 <- expm1(b)
  e2 <- estimate(b, function(p) expm1(p))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  # pow, sqrt
  e1 <- sqrt(e)
  e2 <- estimate(e, function(p) sqrt(p))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  # trigonometric
  for (f in c(cos, sin, tan, acos, asin, atan)) {
      e1 <- f(b)
      e2 <- estimate(b, function(p) f(p))
      testthat::expect_equivalent(IC(e1), IC(e2))
      testthat::expect_equivalent(coef(e1), coef(e2))
  }
  e1 <- sin(asin(b))
  testthat::expect_equivalent(IC(e1), IC(b))
  testthat::expect_equivalent(coef(e1), coef(b))

  # hyberbolic
  for (f in c(cosh, sinh, tanh)) {
      e1 <- f(b)
      e2 <- estimate(b, function(p) f(p))
      testthat::expect_equivalent(IC(e1), IC(e2))
      testthat::expect_equivalent(coef(e1), coef(e2))
  }
  e1 <- acosh(cosh(b))
  testthat::expect_equivalent(IC(e1), IC(b))
  testthat::expect_equivalent(coef(e1), coef(b))
  e1 <- asinh(sinh(b))
  testthat::expect_equivalent(IC(e1), IC(b))
  testthat::expect_equivalent(coef(e1), coef(b))
  e1 <- atanh(tanh(b))
  testthat::expect_equivalent(IC(e1), IC(b))
  testthat::expect_equivalent(coef(e1), coef(b))

})

test_that("custom functions", {
  e1 <- logit(expit(a))
  testthat::expect_equivalent(IC(e1), IC(a))
  testthat::expect_equivalent(coef(e1), coef(a))

  e1 <- 2 * log(a["a2"]) * sqrt(b["b2"])
  e2 <- estimate(merge(a,b),
                 function(p) 2 * log(p["a2"]) * sqrt(p["b2"]))
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  e1 <- 1 / cos(a) - sqrt(a^2 + 1) + 1
  e2 <- estimate(a,
                 function(p) 1 / cos(p) - sqrt(p^2 + 1) + 1)
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))

  e1 <- odds(b)
  e2 <- estimate(b, odds)
  testthat::expect_equivalent(IC(e1), IC(e2))
  testthat::expect_equivalent(coef(e1), coef(e2))
})
