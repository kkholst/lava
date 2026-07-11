## Tests for soft deprecation of null/contrast/type/var.adj in estimate.default

context("estimate.default deprecation warnings")

make_glm <- function() {
  set.seed(1)
  d <- data.frame(
    y = rbinom(100, 1, 0.5),
    x = rnorm(100),
    z = rnorm(100)
  )
  glm(y ~ x + z, data = d, family = binomial())
}

test_that("estimate(..., null=) emits deprecation warning", {
  g <- make_glm()
  expect_warning(estimate(g, null = 0), class = "deprecatedWarning")
})

test_that("estimate(..., back.transform=) emits deprecation warning", {
  g <- make_glm()
  expect_warning(
    ss <- estimate(g, back.transform = exp),
    class = "deprecatedWarning"
  )

  expect_equal(exp(coef(g)), coef(ss))
})

test_that("estimate(..., level=) emits deprecation warning", {
  g <- make_glm()
  s0 <- estimate(g)
  expect_warning(
    ss <- estimate(g, level = 0.5),
    class = "deprecatedWarning"
  )
  expect_true(all(confint(s0) != confint(ss)))
})


test_that("estimate(..., null=0) warns even when null=0 (Q7)", {
  g <- make_glm()
  expect_warning(estimate(g, null = 0), class = "deprecatedWarning")
})

test_that("estimate(..., contrast=) emits deprecation warning", {
  g <- make_glm()
  B <- rbind(c(1, 0, 0), c(0, 1, 0))
  expect_warning(estimate(g, contrast = B), class = "deprecatedWarning")
})

test_that("estimate(..., type=) emits deprecation warning", {
  g <- make_glm()
  expect_warning(estimate(g, type = "df"), class = "deprecatedWarning")
})

test_that("estimate(..., var.adj=) emits deprecation warning", {
  g <- make_glm()
  expect_warning(
    estimate(g, type = "hc3", var.adj = 0.5),
    class = "deprecatedWarning"
  )
})

test_that("multiple deprecated args produce exactly one deprecation warning", {
  g <- make_glm()
  ws <- capture_warnings(estimate(g, null = 1, type = "df"))
  dep_msgs <- grep("deprecated", ws, value = TRUE)
  expect_equal(length(dep_msgs), 1L)
})

test_that("internally-derived contrast from f does NOT warn (Q8)", {
  g <- make_glm()
  ws <- capture_warnings(estimate(g, "x"))
  expect_equal(length(grep("deprecated", ws)), 0L)
  ws <- capture_warnings(estimate(g, 2:3))
  expect_equal(length(grep("deprecated", ws)), 0L)
  ws <- capture_warnings(estimate(g, rbind(c(1, 0, 0), c(0, 1, 0))))
  expect_equal(length(grep("deprecated", ws)), 0L)
})

test_that("plain estimate() emits no deprecation warning (Q4)", {
  g <- make_glm()
  ws <- capture_warnings(estimate(g))
  expect_equal(length(grep("deprecated", ws)), 0L)
})
