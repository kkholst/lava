set.seed(1)
center_ic <- function(n, p = 1) {
  x <- matrix(rnorm(n * p), n, p)
  scale(x, center = TRUE, scale = FALSE)
}
a1 <- estimate(coef = 1, IC = center_ic(10), id = 1:10, labels = "a1")
a2 <- estimate(coef = 2, IC = center_ic(10), id = 1:10, labels = "a2")
a3 <- estimate(coef = 3, IC = center_ic(10), id = 1:10, labels = "a3")

test_that("c.summary.estimate concatenates coefficient matrices", {
  # test basic functionality without any extra variables
  s1 <- summary(a1)
  s2 <- summary(a2)
  s3 <- summary(a3)
  cc <- c(s1, s2, s3)

  # Result is a summary.estimate object
  expect_s3_class(cc, "summary.estimate")

  # coefmat is the row-bind of the inputs (preserving rows and columns)
  expect_equal(rownames(cc$coefmat), c("a1", "a2", "a3"))
  expect_equal(colnames(cc$coefmat), colnames(s1$coefmat))
  expect_equivalent(cc$coefmat, rbind(s1$coefmat, s2$coefmat, s3$coefmat))

  # coef.summary.estimate returns the combined coef vector
  expect_equal(coef(cc), cc$coefmat[, 1])
  expect_equal(parameter(cc), cc$coefmat)

  # original summaries are retained in the objects element
  expect_length(cc$objects, 3L)
  expect_equal(cc$objects[[2]]$coefmat, s2$coefmat)

  # check variance covariance matrix
  expect_equal(dim(vcov(cc)), c(3L, 3L))
})

test_that("c.summary.estimate with single argument returns input unchanged", {
  s1 <- summary(a1)
  expect_identical(c(s1), s1)
})

test_that("c.summary.estimate carries a custom print method", {
  cc <- c(summary(a1), summary(a2))
  expect_type(cc$print, "closure")
  out <- capture.output(print(cc))
  expect_true(any(grepl("Concatenated summary.estimate objects", out)))
  expect_true(any(grepl("a1", out)))
  expect_true(any(grepl("a2", out)))
})

test_that("c.summary.estimate with extra variables", {
  s1 <- summary(a1)
  ss <- c(s1, niter = 1) # single numerical argument
  expect_equal(attributes(ss)$extra, c(niter = 1))

  ss <- c(s1, 1) # single unnamed variable
  expect_equal(attributes(ss)$extra, 1)

  expect_equal(attributes(c(ss, a = 1))$extra, c(1, a = 1))

  # with single named vector
  ss <- c(s1, c(niter = 1, cc = 2)) # single numerical argument
  expect_equal(attributes(ss)$extra, c(niter = 1, cc = 2))

  ss <- c(s1, niter = 1, cc = 2) # multiple numerical arguments
  expect_equal(attributes(ss)$extra, c(niter = 1, cc = 2))

  # two summary.estimate objects, of which one has extra attributes
  ss <- c(s1, ss)
  expect_equal(attributes(ss)$extra, c(niter = 1, cc = 2))

  # two summary.estimate objects, with both having extra attributes
  ss <- c(c(s1, aa = 1), c(s1, bb = 2))
  expect_equal(attributes(ss)$extra, c(aa = 1, bb = 2))

  # adding additional extra to existing extra attribute
  ss <- c(ss, c(cc = 3))
  expect_equal(attributes(ss)$extra, c(aa = 1, bb = 2, cc = 3))
})

test_that("summary.estimate class tests", {
  # tests for future regression when removing estimate class "inheritance"
  # from summary.estimate object class definition
  ss <- summary(a1)
  expect_true(
    all(c(
      inherits(ss, "summary.estimate"),
      inherits(ss, "estimate")
    ))
  )

  # same for concatenation
  ss <- c(summary(a1), aa = 1)
  expect_true(
    all(c(
      inherits(ss, "summary.estimate"),
      inherits(ss, "estimate")
    ))
  )
})
