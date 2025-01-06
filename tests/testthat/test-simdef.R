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
