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
  res <- sim(onerun, 100)
  s <- summary(res)
  expect_true(ncol(s) == 2L)
  expect_equivalent(colnames(s), c("a", "b"))
  expect_equivalent(s["SE",], colMeans(res[, c("Std.Err.a", "Std.Err.b")]))
  expect_equivalent(s["SD",], c(sd(res[,"a"]), sd(res[,"b"])))
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
