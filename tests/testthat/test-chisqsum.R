context("Weighted sum of chi-squared (chisqsum)")

## rchisqsum / pchisqsum are internal -> accessed via lava:::

test_that("rchisqsum returns numeric vector of correct length and non-negative", {
    set.seed(1)
    y <- lava:::rchisqsum(100, lambda = c(1, 2, 0.5))
    testthat::expect_type(y, "double")
    testthat::expect_length(y, 100)
    testthat::expect_true(all(y >= 0))
})

test_that("rchisqsum mean and variance match theoretical values", {
    testthat::skip_on_cran()
    set.seed(42)
    lambda <- c(0.5, 1, 2)
    y <- lava:::rchisqsum(2e5, lambda = lambda)
    ## E[sum lambda_i Z_i^2] = sum(lambda)
    testthat::expect_lt(abs(mean(y) - sum(lambda)), 0.05)
    ## Var = 2 * sum(lambda^2)
    testthat::expect_lt(abs(var(y) - 2 * sum(lambda^2)), 0.5)
})

test_that("rchisqsum reduces to chi-squared(df) when lambda = rep(1, df)", {
    testthat::skip_on_cran()
    set.seed(7)
    df <- 5
    y <- lava:::rchisqsum(2e5, lambda = rep(1, df))
    testthat::expect_lt(abs(mean(y) - df), 0.05)
    testthat::expect_lt(abs(var(y) - 2 * df), 0.5)
})

test_that("rchisqsum with single lambda=1 matches chi-squared(1)", {
    testthat::skip_on_cran()
    set.seed(11)
    y <- lava:::rchisqsum(2e5, lambda = 1)
    testthat::expect_lt(abs(mean(y) - 1), 0.03)
    testthat::expect_lt(abs(var(y) - 2), 0.1)
})

test_that("pchisqsum returns probability in [0,1]", {
    set.seed(1)
    p <- lava:::pchisqsum(1, lambda = 1, B = 5e3, seed = 1)
    testthat::expect_true(p >= 0 && p <= 1)
})

test_that("pchisqsum is monotone non-decreasing in x", {
    testthat::skip_on_cran()
    lambda <- c(1, 2)
    xs <- c(0.5, 1, 2, 5, 10)
    ps <- vapply(xs,
                 function(xx) lava:::pchisqsum(xx, lambda = lambda,
                                                B = 2e4, seed = 1),
                 numeric(1))
    testthat::expect_true(all(diff(ps) >= 0))
})

test_that("pchisqsum boundary behaviour at extremes", {
    testthat::skip_on_cran()
    p_lo <- lava:::pchisqsum(-1, lambda = 1, B = 1e4, seed = 1)
    p_hi <- lava:::pchisqsum(1e6, lambda = 1, B = 1e4, seed = 1)
    testthat::expect_equal(p_lo, 0)
    testthat::expect_equal(p_hi, 1)
})

test_that("pchisqsum is reproducible with seed argument", {
    p1 <- lava:::pchisqsum(2, lambda = c(1, 0.5), B = 5e3, seed = 99)
    p2 <- lava:::pchisqsum(2, lambda = c(1, 0.5), B = 5e3, seed = 99)
    testthat::expect_identical(p1, p2)
    p3 <- lava:::pchisqsum(2, lambda = c(1, 0.5), B = 5e3, seed = 100)
    testthat::expect_false(isTRUE(all.equal(p1, p3)))
})

test_that("pchisqsum approximates pchisq when lambda = 1", {
    testthat::skip_on_cran()
    xs <- c(0.5, 1, 2, 4)
    for (x in xs) {
        p_mc <- lava:::pchisqsum(x, lambda = 1, B = 2e5, seed = 13)
        p_exact <- pchisq(x, df = 1)
        testthat::expect_lt(abs(p_mc - p_exact), 0.01)
    }
})

## ---------------------------------------------------------------------------
## Regression: pchisqsum vectorized over x
## ---------------------------------------------------------------------------
test_that("pchisqsum returns one probability per element of vector x", {
    xs <- c(0.5, 1, 2, 4)
    p_vec <- lava:::pchisqsum(xs, lambda = c(1, 0.5), B = 5e3, seed = 7)
    testthat::expect_length(p_vec, length(xs))
    testthat::expect_true(all(p_vec >= 0 & p_vec <= 1))
    testthat::expect_true(all(diff(p_vec) >= 0))
    ## scalar call agrees with vectorized call (same seed, same B)
    p1 <- lava:::pchisqsum(xs[2], lambda = c(1, 0.5), B = 5e3, seed = 7)
    testthat::expect_equal(p_vec[2], p1)
})
