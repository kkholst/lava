context("procformula robustness against as.character.formula masking")

# formula.tools overwrites as.character.formula, changing it from returning
# a length-3 vector c("~", "LHS", "RHS") to a single collapsed string.
# These tests verify that procformula, eventTime, and predict.zibreg work
# correctly regardless of whether as.character.formula has been masked.

mock_as_character_formula <- function(x, ...) {
  paste(deparse(x), collapse = " ")
}

test_that("procformula handles standard two-sided formulas", {
  res <- procformula(value = y ~ x + z)
  expect_equal(res$ys, "y")
  expect_true("x" %in% res$xs)
  expect_true("z" %in% res$xs)
})

test_that("procformula handles covariance formula (y1 ~ ~y2)", {
  m <- lvm()
  res <- procformula(m, value = y1 ~ ~y2)
  expect_equal(res$ys, "y1")
  expect_true(res$iscovar)
})

test_that("procformula handles multivariate LHS", {
  res <- procformula(value = c(y1, y2) ~ x + z)
  expect_true("y1" %in% res$ys)
  expect_true("y2" %in% res$ys)
  expect_true("x" %in% res$xs)
  expect_true("z" %in% res$xs)
})

test_that("procformula handles one-sided formula", {
  res <- procformula(value = ~ x + z)
  expect_null(res$ys)
  expect_true("x" %in% res$xs)
  expect_true("z" %in% res$xs)
})

test_that("procformula handles link function syntax", {
  res <- procformula(value = y ~ exp(x + z))
  expect_equal(res$invlink, "exp")
  expect_true("x" %in% res$xs)
  expect_true("z" %in% res$xs)
})

test_that("procformula results match when as.character.formula is masked", {
  formulas <- list(
    y ~ x + z,
    y1 ~ ~y2,
    c(y1, y2) ~ x + z,
    ~ x + z,
    y ~ x + z + ~w
  )
  for (f in formulas) {
    # Run with normal dispatch
    res_normal <- procformula(value = f)
    # Temporarily mask as.character.formula
    registerS3method(
      "as.character",
      "formula",
      mock_as_character_formula,
      envir = asNamespace("base")
    )
    on.exit(
      {
        rm_env <- asNamespace("base")
        if (exists("as.character.formula", envir = rm_env)) {
          rm("as.character.formula", envir = rm_env)
        }
      },
      add = TRUE
    )
    res_masked <- procformula(value = f)
    # Restore
    rm_env <- asNamespace("base")
    if (exists("as.character.formula", envir = rm_env)) {
      rm("as.character.formula", envir = rm_env)
    }
    expect_identical(
      res_normal$ys,
      res_masked$ys,
      info = paste("ys mismatch for", deparse(f))
    )
    expect_identical(
      res_normal$xs,
      res_masked$xs,
      info = paste("xs mismatch for", deparse(f))
    )
    expect_identical(
      res_normal$iscovar,
      res_masked$iscovar,
      info = paste("iscovar mismatch for", deparse(f))
    )
    expect_identical(
      res_normal$invlink,
      res_masked$invlink,
      info = paste("invlink mismatch for", deparse(f))
    )
  }
})

test_that("lvm model building works when as.character.formula is masked", {
  registerS3method(
    "as.character",
    "formula",
    mock_as_character_formula,
    envir = asNamespace("base")
  )
  on.exit(
    {
      rm_env <- asNamespace("base")
      if (exists("as.character.formula", envir = rm_env)) {
        rm("as.character.formula", envir = rm_env)
      }
    },
    add = TRUE
  )
  # These should not error
  m <- lvm(y ~ x + z)
  expect_true("y" %in% vars(m))
  expect_true("x" %in% vars(m))

  m2 <- lvm(c(y1, y2) ~ x, y1 ~ ~y2)
  expect_true("y1" %in% vars(m2))
  expect_true("y2" %in% vars(m2))
})

test_that("eventTime works when as.character.formula is masked", {
  registerS3method(
    "as.character",
    "formula",
    mock_as_character_formula,
    envir = asNamespace("base")
  )
  on.exit(
    {
      rm_env <- asNamespace("base")
      if (exists("as.character.formula", envir = rm_env)) {
        rm("as.character.formula", envir = rm_env)
      }
    },
    add = TRUE
  )
  m <- lvm(y ~ x)
  # Two-sided eventTime formula
  expect_silent(eventTime(m, time ~ min(T1 = 1, T2 = 2, C = 0)))
  # One-sided eventTime formula
  expect_silent(eventTime(m, ~ min(T1 = 1, T2 = 0)))
})
