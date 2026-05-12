## Equivalence tests: estimate(g, X=...) should produce the same coefmat/vcov
## as summary(estimate(g), X=...) during the soft-deprecation window.

context("estimate / summary equivalence for migrated args")

make_glm <- function() {
  set.seed(1)
  d <- data.frame(
    y = rbinom(100, 1, 0.5),
    x = rnorm(100),
    z = rnorm(100)
  )
  glm(y ~ x + z, data = d, family = binomial())
}

## Helper: get coefmat from either deprecated or new path
deprecated_estimate <- function(...) {
  suppressWarnings(estimate(...))
}

test_that("null= gives identical coefmat via estimate() and summary()", {
  g <- make_glm()
  e1 <- deprecated_estimate(g, null = 1)
  e2 <- summary(estimate(g), null = 1)
  expect_equal(e1$coefmat, e2$coefmat)
  expect_equal(unname(e1$vcov), unname(e2$vcov))
})

test_that("contrast= gives identical coefmat via estimate() and summary()", {
  g <- make_glm()
  B <- rbind(c(1, 0, 0), c(0, 1, -1))
  e1 <- deprecated_estimate(g, contrast = B)
  e2 <- summary(estimate(g), contrast = B)
  expect_equal(e1$coefmat, e2$coefmat)
  expect_equal(unname(e1$vcov), unname(e2$vcov))
})

test_that("null + contrast gives identical results via estimate() and summary()", {
  g <- make_glm()
  B <- rbind(c(1, 0, 0), c(0, 1, -1))
  e1 <- deprecated_estimate(g, contrast = B, null = c(0, 0))
  e2 <- summary(estimate(g), contrast = B, null = c(0, 0))
  expect_equal(e1$coefmat, e2$coefmat)
  expect_equal(unname(e1$vcov), unname(e2$vcov))
})

test_that("type='df' gives identical coefmat via estimate() and summary()", {
  g <- make_glm()
  e1 <- deprecated_estimate(g, type = "df")
  e2 <- summary(estimate(g), type = "df")
  expect_equal(e1$coefmat, e2$coefmat)
  expect_equal(unname(e1$vcov), unname(e2$vcov))
})

test_that("type='df1' gives identical coefmat via estimate() and summary()", {
  g <- make_glm()
  e1 <- deprecated_estimate(g, type = "df1")
  e2 <- summary(estimate(g), type = "df1")
  expect_equal(e1$coefmat, e2$coefmat)
})

test_that("type='hc3' + var.adj gives identical results via estimate() and summary()", {
  g <- make_glm()
  e1 <- deprecated_estimate(g, type = "hc3", var.adj = 0.5)
  e2 <- summary(estimate(g), type = "hc3", var.adj = 0.5)
  expect_equal(e1$coefmat, e2$coefmat)
  expect_equal(unname(e1$vcov), unname(e2$vcov))
})

test_that("type + null gives identical results via estimate() and summary()", {
  g <- make_glm()
  e1 <- deprecated_estimate(g, type = "df", null = 0.2)
  e2 <- summary(estimate(g), type = "df", null = 0.2)
  expect_equal(e1$coefmat, e2$coefmat)
})

test_that("df + null gives identical results via estimate() and summary()", {
  g <- make_glm()
  ## Per Q11: df is stored on the estimate object and reused by summary().
  e1 <- deprecated_estimate(g, df = 10, null = 1)
  e2 <- summary(estimate(g, df = 10), null = 1)
  expect_equal(e1$coefmat, e2$coefmat)
})

test_that("df is stored on the estimate object", {
  g <- make_glm()
  e <- estimate(g, df = 7)
  expect_equal(e$df, 7)
  e0 <- estimate(g)
  expect_null(e0$df)
})
