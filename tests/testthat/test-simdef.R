library(future)
context("Generic simulation framework")

f <- function(iter=list(), ...) {
  if (!is.list(iter) || is.null(iter$i)) return(0)
    return(iter$i)
}
R <- list(
  list(i = 2),
  list(i = 4)
)

test_that("sim.default, list in put", {
  res <- sim(f, R)
  testthat::expect_true(length(R) == nrow(res))
  testthat::expect_identical(
    as.vector(res),
    unname(unlist(R))
  )
})

test_that("sim.default with estimate objects", {
  onerun <- function(...) estimate(coef=runif(2),
                                   vcov=diag(runif(2)),
                                   labels=c("a","b"))
  res <- sim(onerun, 10)
  s <- summary(res)
  expect_true(ncol(s) == 2L)
  expect_equivalent(colnames(s), c("a", "b"))
  expect_equivalent(s["SE",], colMeans(res[, c("a.Std.Err", "b.Std.Err")]))
  expect_equivalent(s["SD",], c(sd(res[,"a"]), sd(res[,"b"])))
})

test_that("sim.default with summary.estimate objects", {
  onerun <- function(...) {
    estimate(coef = runif(2), vcov = diag(runif(2)), labels = c("a","b")) |>
      summary()
  }
  res <- sim(onerun, 5)
  s <- summary(res)
  expect_true(ncol(s) == 2L)
  expect_equivalent(colnames(s), c("a", "b"))
  expect_equivalent(s["SE",], colMeans(res[, c("a.Std.Err", "b.Std.Err")]))
  expect_equivalent(s["SD",], c(sd(res[,"a"]), sd(res[,"b"])))

  # with concatenated summary.estimate objects
  onerun <- function(...) {
    s1 <- estimate(
      coef = runif(2), vcov = diag(runif(2)), labels = c("a","b")
    ) |> summary()
    s2 <- estimate(
      coef = runif(2), vcov = diag(runif(2)), labels = c("c","d")
    ) |> summary(transform = exp) # blanks Std.Err column in coefmat
    c(s1, s2)
  }
  res <- sim(onerun, 5)
  expect_true(all(is.na(res[, c("c.Std.Err", "d.Std.Err")])))

  s <- summary(res)
  expect_true(ncol(s) == 4L)
  expect_equivalent(colnames(s), c("a", "b", "c", "d"))

  # with concatenated summary.estimate objects + extra
  onerun <- function(...) {
    s1 <- estimate(
      coef = runif(2), vcov = diag(runif(2)), labels = c("a","b")
    ) |> summary()
    c(s1, niter = rpois(1, 10), converged = 1)
  }
  res <- sim(onerun, 5)
  expect_true(all(c("niter", "converged") %in% colnames(res)))
  expect_true(all(res[, "converged"] == 1))
})

test_that("sim.default exports seed sequences as attribute", {
  foo <- function() runif(1)
  if (requireNamespace("future",quietly=TRUE)) future::plan("sequential")
  result <- sim(foo, R = 5, future.seed = 42L)
  seeds <- attr(result, "seeds")

  expect_true(is.list(seeds))
  expect_equal(length(seeds), 5L)
  expect_true(is.integer(seeds[[1]]))
})

test_that("sim.default exported seeds reproduce results (sequential)", {
  foo <- function() runif(1)
  if (requireNamespace("future",quietly=TRUE)) future::plan("sequential")
  result <- sim(foo, R = 5, future.seed = 42L)
  seeds <- attr(result, "seeds")

  old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))

  for (i in seq_len(5)) {
    assign(".Random.seed", seeds[[i]], envir = .GlobalEnv)
    expect_equal(foo(), as.numeric(result[i, 1]))
  }
})

test_that("sim.default exported seeds reproduce results (mc.cores = 1)", {
  skip_on_os("windows")
  foo <- function() runif(1)
  result <- sim(foo, R = 5, mc.cores = 1L)
  seeds <- attr(result, "seeds")

  old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))

  for (i in seq_len(5)) {
    assign(".Random.seed", seeds[[i]], envir = .GlobalEnv)
    expect_equal(foo(), as.numeric(result[i, 1]))
  }
})

test_that("sim.default subsets", {
  onerun <- function(...) estimate(coef=runif(2),
                                   vcov=diag(runif(2)),
                                   labels=c("a","b"))
  res <- sim(onerun, 100, estimate.index = 1)

  expect_true(ncol(res) == 2L)
  expect_true(nrow(res) == 100L)

  expect_true(nrow(
    rbind(res, res)
  ) == 200L)

  expect_true(nrow(
    rbind(res[1:20, 1], res[1:50, 1])
  ) == 70L)

  expect_true(ncol(
    cbind(res[1:50, 1], res[1:50, 2])
  ) == 2L)

  expect_true(ncol(
    cbind(res[1:50, 1:2], res[1:50, 1:2])
  ) == 4L)

})

test_that("sim.default with estimate.extra appends extras", {
  onerun <- function(...) {
    e <- estimate(coef = runif(2), vcov = diag(runif(2)),
                  labels = c("a", "b"))
    c(e, converged = 1, niter = sample(5:20, 1))
  }
  res <- sim(onerun, 20)
  pi <- attr(res, "par.index")

  # Correct dimensions: 2 estimates + 2 SE + 2 extras = 6 columns
  expect_equal(ncol(res), 6L)
  expect_equal(nrow(res), 20L)

  # par.index tracks all components

  expect_equal(pi$estimate, 1:2)
  expect_equal(pi$se, 3:4)
  expect_equal(pi$extra, 5:6)

  # Extra columns have correct names
  expect_true("converged" %in% colnames(res))
  expect_true("niter" %in% colnames(res))

  # All converged values are 1
  expect_true(all(res[, "converged"] == 1))
  # niter values are in expected range
  expect_true(all(res[, "niter"] >= 5 & res[, "niter"] <= 20))
})

test_that("summary.sim ignores extra columns", {
  onerun <- function(...) {
    e <- estimate(coef = runif(2), vcov = diag(runif(2)),
                  labels = c("a", "b"))
    c(e, flag = 1)
  }
  res <- sim(onerun, 30)
  s <- summary(res)

  # Summary only covers the estimate columns, not extras
  expect_equal(ncol(s), 2L)
  expect_equal(colnames(s), c("a", "b"))
})

test_that("summary.sim gives NA SE for extras when estimate is explicit", {
  onerun <- function(...) {
    e <- estimate(coef = runif(2), vcov = diag(runif(2)),
                  labels = c("a", "b"))
    c(e, conv = 1)
  }
  res <- sim(onerun, 20)
  # Include extra column in explicit estimate argument
  s <- summary(res, estimate = c(1, 2, 5))

  expect_equal(ncol(s), 3L)
  # SE is NA for the extra column
  expect_true(is.na(s["SE", "conv"]))
  expect_true(is.na(s["SE/SD", "conv"]))
  # SE is computed for the estimate columns
  expect_false(is.na(s["SE", "a"]))
  expect_false(is.na(s["SE", "b"]))
})
