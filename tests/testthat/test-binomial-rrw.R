context("Binomial risk-difference / relative-risk models (binomial.rrw)")

## ---------------------------------------------------------------------------
## Identical()
## ---------------------------------------------------------------------------
test_that("Identical: numerical equality with default tolerance", {
  testthat::expect_true(Identical(1, 1))
  testthat::expect_true(Identical(1 + 1e-12, 1))
  testthat::expect_false(Identical(1 + 1e-3, 1))
  ## vector input
  res <- Identical(c(1, 2, 1 + 1e-12), 1)
  testthat::expect_equal(res, c(TRUE, FALSE, TRUE))
  ## custom tolerance
  testthat::expect_true(Identical(1.001, 1, tolerance = 1e-2))
})

## ---------------------------------------------------------------------------
## RD_OP: round-trip check
## ---------------------------------------------------------------------------
test_that("RD_OP recovers (p0, p1) from (rd, op)", {
  RD_OP <- lava:::RD_OP
  p0 <- c(0.2, 0.4, 0.6)
  p1 <- c(0.3, 0.5, 0.8)
  rd <- p1 - p0
  ## op is the odds *product*: (p0/(1-p0)) * (p1/(1-p1))
  op <- (p0 / (1 - p0)) * (p1 / (1 - p1))
  out <- RD_OP(rd, op)
  testthat::expect_equal(dim(out), c(3L, 2L))
  testthat::expect_equal(out[, 1], p0, tolerance = 1e-8)
  testthat::expect_equal(out[, 2], p1, tolerance = 1e-8)
})

test_that("RD_OP handles op == 1 (no association) branch", {
  RD_OP <- lava:::RD_OP
  rd <- c(0.1, -0.2)
  op <- c(1, 1)
  out <- RD_OP(rd, op)
  ## When op == 1: p0 = 0.5*(1 - rd), p1 = p0 + rd
  testthat::expect_equal(out[, 1], 0.5 * (1 - rd), tolerance = 1e-10)
  testthat::expect_equal(out[, 2], 0.5 * (1 - rd) + rd, tolerance = 1e-10)
  testthat::expect_true(all(is.finite(out)))
})

## ---------------------------------------------------------------------------
## RR_OP: round-trip check
## ---------------------------------------------------------------------------
test_that("RR_OP recovers (p0, p1) from (rr, op)", {
  RR_OP <- lava:::RR_OP
  p0 <- c(0.2, 0.3, 0.5)
  p1 <- c(0.3, 0.45, 0.6)
  rr <- p1 / p0
  ## op is the odds *product*: (p0/(1-p0)) * (p1/(1-p1))
  op <- (p0 / (1 - p0)) * (p1 / (1 - p1))
  out <- RR_OP(rr, op)
  testthat::expect_equal(dim(out), c(3L, 2L))
  testthat::expect_equal(out[, 1], p0, tolerance = 1e-8)
  testthat::expect_equal(out[, 2], p1, tolerance = 1e-8)
})

test_that("RR_OP handles op == 1 (no association) branch", {
  RR_OP <- lava:::RR_OP
  rr <- c(1.5, 2)
  op <- c(1, 1)
  out <- RR_OP(rr, op)
  ## When op == 1: p0 = 1/(1+rr), p1 = p0 * rr
  testthat::expect_equal(out[, 1], 1 / (1 + rr), tolerance = 1e-10)
  testthat::expect_equal(out[, 2], rr / (1 + rr), tolerance = 1e-10)
  testthat::expect_true(all(is.finite(out)))
})

## ---------------------------------------------------------------------------
## binomial.rd model modifier
## ---------------------------------------------------------------------------
test_that("binomial.rd installs multiple-input simulation node", {
  m <- lvm()
  regression(m) <- z ~ x
  regression(m) <- lp ~ x
  regression(m) <- op ~ x
  m <- binomial.rd(
    m,
    response = "y",
    exposure = "z",
    target.model = "lp",
    nuisance.model = "op"
  )
  mi <- m$attributes$multiple.inputs
  testthat::expect_true(!is.null(mi))
  testthat::expect_true("y" %in% names(mi))
  testthat::expect_equal(mi$y$input, c("z", "lp", "op"))
  testthat::expect_true(grepl("risk-difference", mi$y$type, fixed = TRUE))
  ## exposure should now have a binomial distribution attached
  testthat::expect_false(is.null(distribution(m)[["z"]]))
})

test_that("binomial.rd accepts formula interface", {
  m1 <- lvm()
  regression(m1) <- z ~ x
  regression(m1) <- lp ~ x
  regression(m1) <- op ~ x
  m1 <- binomial.rd(
    m1,
    response = "y",
    exposure = "z",
    target.model = "lp",
    nuisance.model = "op"
  )

  m2 <- lvm()
  regression(m2) <- z ~ x
  regression(m2) <- lp ~ x
  regression(m2) <- op ~ x
  m2 <- binomial.rd(m2, y ~ z | lp | op)

  testthat::expect_equal(
    m1$attributes$multiple.inputs$y$input,
    m2$attributes$multiple.inputs$y$input
  )
  testthat::expect_equal(
    m1$attributes$multiple.inputs$y$type,
    m2$attributes$multiple.inputs$y$type
  )
})

## ---------------------------------------------------------------------------
## binomial.rr model modifier
## ---------------------------------------------------------------------------
test_that("binomial.rr installs multiple-input simulation node", {
  m <- lvm()
  regression(m) <- z ~ x
  regression(m) <- lp ~ x
  regression(m) <- op ~ x
  m <- binomial.rr(
    m,
    response = "y",
    exposure = "z",
    target.model = "lp",
    nuisance.model = "op"
  )
  mi <- m$attributes$multiple.inputs
  testthat::expect_true("y" %in% names(mi))
  testthat::expect_equal(mi$y$input, c("z", "lp", "op"))
  testthat::expect_true(grepl("relative-risk", mi$y$type, fixed = TRUE))
})

## ---------------------------------------------------------------------------
## End-to-end simulation: numerical correctness
## ---------------------------------------------------------------------------
test_that("binomial.rd simulation produces target risk difference", {
  testthat::skip_on_cran()
  ## Build model where lp is constant (intercept-only) so the implied RD
  ## is tanh(lp_const).  We set lp via constrain().
  m <- lvm()
  regression(m) <- z ~ x
  distribution(m, ~x) <- normal.lvm()
  distribution(m, ~lp) <- normal.lvm(mean = 0, sd = 0) ## degenerate -> 0
  distribution(m, ~op) <- normal.lvm(mean = 0, sd = 0) ## op = exp(0) = 1
  intercept(m, ~lp) <- 0.4 ## constant linear predictor for RD
  intercept(m, ~op) <- 0 ## odds product = 1
  m <- binomial.rd(
    m,
    response = "y",
    exposure = "z",
    target.model = "lp",
    nuisance.model = "op"
  )
  set.seed(123)
  d <- sim(m, n = 8000)
  testthat::expect_true(all(d$y %in% c(0, 1)))
  rd_emp <- mean(d$y[d$z == 1]) - mean(d$y[d$z == 0])
  target <- tanh(0.4)
  ## SE roughly sqrt(p*(1-p)/n_per_group) ~ 0.01; allow generous 3x
  testthat::expect_lt(abs(rd_emp - target), 0.04)
})

test_that("binomial.rr simulation produces target relative risk", {
  testthat::skip_on_cran()
  m <- lvm()
  regression(m) <- z ~ x
  distribution(m, ~x) <- normal.lvm()
  distribution(m, ~lp) <- normal.lvm(mean = 0, sd = 0)
  distribution(m, ~op) <- normal.lvm(mean = 0, sd = 0)
  intercept(m, ~lp) <- log(1.5) ## constant log-RR
  intercept(m, ~op) <- 0
  m <- binomial.rr(
    m,
    response = "y",
    exposure = "z",
    target.model = "lp",
    nuisance.model = "op"
  )
  set.seed(456)
  d <- sim(m, n = 8000)
  testthat::expect_true(all(d$y %in% c(0, 1)))
  p1 <- mean(d$y[d$z == 1])
  p0 <- mean(d$y[d$z == 0])
  rr_emp <- p1 / p0
  testthat::expect_lt(abs(log(rr_emp) - log(1.5)), 0.10)
})

## ---------------------------------------------------------------------------
## Bonus: riskcomp (lives in R/assoc.R but logically related)
## ---------------------------------------------------------------------------
test_that("riskcomp computes risk difference / risk ratio for two probabilities", {
  p <- c(0.4, 0.6)
  ## Default op="/", type=1: val[1]/val[2]
  testthat::expect_equal(riskcomp(p[1], p[2]), p[1] / p[2])
  ## Explicit Diff (op="-")
  testthat::expect_equal(Diff(p[1], p[2]), p[1] - p[2])
  ## Explicit Ratio (op="/")
  testthat::expect_equal(Ratio(p[1], p[2]), p[1] / p[2])
  ## type=2 swaps operand order
  testthat::expect_equal(riskcomp(p[1], p[2], type = 2), p[2] / p[1])
})

test_that("riskcomp with scale=odds gives odds ratio", {
  p <- c(0.4, 0.6)
  o <- p / (1 - p)
  testthat::expect_equal(riskcomp(p[1], p[2], scale = odds), o[1] / o[2])
  ## type=2 swaps order
  testthat::expect_equal(
    riskcomp(p[1], p[2], scale = odds, type = 2),
    o[2] / o[1]
  )
})

## ---------------------------------------------------------------------------
## Regression: simulate_binomial_* errors clearly on out-of-range probabilities
## ---------------------------------------------------------------------------
test_that("simulate_binomial_rd errors when implied probability is out of [0,1]", {
  sim_rd <- lava:::simulate_binomial_rd
  ## NaN linear predictor propagates to NaN probabilities -> guard triggers
  data <- data.frame(z = c(0, 1), lp = c(NaN, NaN), op = c(0, 0))
  testthat::expect_error(
    sim_rd(NULL, data, inputs = c("z", "lp", "op")),
    regexp = "outside"
  )
})

test_that("simulate_binomial_rr errors when implied probability is out of [0,1]", {
  sim_rr <- lava:::simulate_binomial_rr
  data <- data.frame(z = c(0, 1), lp = c(NaN, NaN), op = c(0, 0))
  testthat::expect_error(
    sim_rr(NULL, data, inputs = c("z", "lp", "op")),
    regexp = "outside"
  )
})
