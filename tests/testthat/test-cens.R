context("Censored observations")

library("survival")

test_that("right censoring", {
  set.seed(1)
  n <- 500
  x <- rnorm(n)
  y0 <- rnorm(n) + x
  rcens <- rexp(n, 1)
  y <- pmin(y0, rcens)
  rstatus <- y0 < rcens
  d <- data.frame(
    y = y,
    x = x,
    status = rstatus
  )
  d$ys <- with(d, Surv(y, status, type = "right"))

  ms <- lvm(ys ~ x)
  e <- estimate(ms, data = d)
  es <- survreg(ys ~ x, data = d, dist = "gaussian")

  expect_true(abs(logLik(e) - logLik(es)) < 1e-3)
  expect_true(mean((coef(e)[1:2] - coef(es))^2) < 1e-3)
})

test_that("left censoring", {
  set.seed(1)
  n <- 500
  x <- rnorm(n)
  y0 <- rnorm(n) + x
  lcens <- -rexp(n, 1)
  y <- pmax(y0, lcens)
  lstatus <- y0 > lcens
  d <- data.frame(
    y = y,
    x = x,
    status = lstatus
  )
  d$ys <- with(d, Surv(y, status, type = "left"))

  ms <- lvm(ys ~ x)
  e <- estimate(ms, data = d)
  es <- survreg(ys ~ x, data = d, dist = "gaussian")

  expect_true(abs(logLik(e) - logLik(es)) < 1e-3)
  expect_true(mean((coef(e)[1:2] - coef(es))^2) < 1e-3)
})

test_that("left and right censoring", {
  set.seed(1)
  n <- 500
  x <- rnorm(n)
  y0 <- rnorm(n) + x
  lcens <- -rexp(n, 2)
  rcens <- rexp(n, 1)
  y <- pmax(y0, lcens)
  y <- pmin(y, rcens)
  lstatus <- y0 > lcens
  rstatus <- y0 < rcens
  d <- data.frame(
    y = y,
    x = x,
    lstatus = lstatus,
    rstatus = rstatus
  )
  yleft <- rep(-Inf, n)
  yleft[d$lstatus] <- d$y[d$lstatus]
  yright <- rep(Inf, n)
  yright[d$rstatus] <- d$y[d$rstatus]
  d$ys <- Surv(yleft, yright, type = "interval2")

  ms <- lvm(ys ~ x)
  e <- estimate(ms, data = d)
  es <- survreg(ys ~ x, data = d, dist = "gaussian")

  expect_true(abs(logLik(e) - logLik(es)) < 1e-3)
  expect_true(mean((coef(e)[1:2] - coef(es))^2) < 1e-3)
})

test_that("interval censoring", {
  set.seed(1)
  n <- 1e4
  x <- rnorm(n) # covariate (no effects.lvmfit)
  true_times <- rweibull(n, shape = 1.5, scale = 50) # true event timestamp
  # Define fixed inspection points every 10 time units
  visits <- seq(0, 100, by = 10)
  # left and right interval (obs: right censored after last visit)
  yleft <- sapply(true_times, function(t) max(visits[visits <= t]))
  yright <- sapply(true_times, function(t) min(visits[visits > t], Inf))

  d <- data.frame(
    x = x,
    yleft = yleft,
    yright = yright
  )
  d$ys <- Surv(yleft, yright, type = "interval2")

  ms <- lvm(ys ~ x)
  e <- estimate(ms, data = d)
  es <- survreg(ys ~ x, data = d, dist = "gaussian")

  expect_true(abs(logLik(e) - logLik(es)) < 1e-3)
  expect_true(mean((coef(e)[1:2] - coef(es))^2) < 1e-3)
})
