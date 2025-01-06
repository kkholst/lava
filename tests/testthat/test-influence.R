context("Influence functions")

test_that("GEE", {
  if (requireNamespace("geepack",quietly=TRUE)) {
  d <- lvm(y ~ x, ~ id) |>
    distribution(~id, uniform.lvm(value=seq(1:20))) |>
    sim(100, seed=1)
  d0 <- d[order(d$id), ]
  g <- geepack::geeglm(y ~ x, data=d0, id=d0$id)
  V <- summary(g)$cov.scaled
  g0 <- glm(y ~ x, data=d)
  V0 <- vcov(estimate(g0, id = d$id))
  testthat::expect_true(sum((V - V0)^2) < 1e-12)
  }
})

test_that("merge, IC, estimate with 'id' argument", {
  d <- data.frame(id=c("a","a","b","b","b","b","c","c","d"),
                id1=c("a","a","b1","b1","b2","b2","c","c","d"),
                y=rnorm(9), x=rnorm(9))
  d$id0 <- as.numeric(as.factor(d$id))

  l <- glm(y ~ x, data=d)
  e1 <- estimate(l, id=d$id1)
  e <- merge(estimate(l), estimate(l), id=list(d$id, d$id1))
  V0 <- vcov(e)
  V1 <- vcov(estimate(e1, id=c(1,2,2,3,4)))

  if (requireNamespace("geepack",quietly=TRUE)) {
    V <- summary(geepack::geeglm(y ~ x, id=d$id0, data=d))$cov.scaled
    testthat::expect_true(sum((V - V0)^2) < 1e-12)
    testthat::expect_true(sum((V - V1)^2) < 1e-12)
  }

  testthat::expect_true(sum((vcov(e1) - vcov(e)[3:4,3:4])^2) < 1e-12)
  testthat::expect_true(sum((V1 - vcov(e)[1:2,1:2])^2) < 1e-12)

  ee <- estimate(e, id=c(1,2,3,4,2,2))
  VV <- vcov(ee)
  testthat::expect_true(sum((VV[1:2,1:2] - V)^2) < 1e-12)
  testthat::expect_true(sum((VV[3:4,3:4] - V)^2) < 1e-12)
})

test_that("negative binomial regression (glm.nb)", {
  set.seed(1)
  n <- 500
  z <- rgamma(n, .5, .5)
  x <- rnorm(n)
  lam <- z * exp(x)
  y <- rpois(n, lam)
  if (requireNamespace("MASS",quietly=TRUE)) {
  m <- MASS::glm.nb(y ~ x)
  testthat::expect_true(abs(lava:::logL.glm(m) - logLik(m)) < 1e-6)
  p <- coef(m)+1
  u1 <- as.vector(numDeriv::jacobian(function(p) lava:::logL.glm(m, p = p), p))
  u2 <- score(m, p = p)
  testthat::expect_true(sum((u1 - u2)^2) < 1e-6)
  p <- coef(m)
  u1 <- as.vector(numDeriv::jacobian(function(p) lava:::logL.glm(m, p = p), p))
  u2 <- score(m, p = p)
  testthat::expect_true(sum((u1 - u2)^2) < 1e-6)
  }
})


test_that("quasipossion", {
  set.seed(1)
  n <- 500
  z <- rgamma(n, .5, .5)
  x <- rnorm(n)
  lam <- z * exp(x)
  y <- rpois(n, lam)
  m1 <- glm(y ~ x, family=poisson)
  m2 <- glm(y ~ x, family = quasipoisson)
  i1 <- IC(m1)
  i2 <- IC(m2)
  testthat::expect_true(sum((i1 - i2)^2) < 1e-6)
})
