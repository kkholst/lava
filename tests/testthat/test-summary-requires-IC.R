## summary.estimate() requires IC=TRUE on the input estimate when applying
## small-sample corrections (type / var.adj). Per Q5.

context("summary.estimate IC requirement")

make_glm <- function() {
  set.seed(1)
  d <- data.frame(
    y = rbinom(100, 1, 0.5),
    x = rnorm(100),
    z = rnorm(100)
  )
  glm(y ~ x + z, data = d, family = binomial())
}

test_that("summary(..., type=) errors when input estimate has IC=FALSE", {
  # variance correction via the type argument requires e$IC to be not null
  g <- make_glm()
  e <- estimate(g, IC = FALSE)
  expect_error(
    summary(e, type = "df"),
    "IC=TRUE",
    fixed = TRUE
  )
})

test_that("summary(..., null=) works when input estimate has IC=FALSE", {
  ## null/contrast do NOT require IC, only type/var.adj do.
  g <- make_glm()
  e <- estimate(g, IC = FALSE)
  expect_silent(summary(e, null = 1))
})
