context("estimate.default id/cluster")

test_that("index.estimate", {
    m <- lvm(c(y1,y2)~x+z, y1~~y2)
    set.seed(1)
    d <- sim(m,20)

    l1 <- lm(y1~x+z, data=d)
    e1 <- estimate(l1)
    testthat::expect_equivalent(index(e1), rownames(d))
    testthat::expect_true(inherits(index(e1), "character"))
    V <- vcov(e1)
    index(e1) <- as.numeric(index(e1))
    expect_equivalent(vcov(e1), V)
    testthat::expect_true(inherits(index(e1), "numeric"))
})


test_that("estimate index order", {
  n <- 20
  d <- data.frame(
    y = rnorm(n),
    a = rbinom(n,1,0.5),
    id = 1:n,
    id2 = n:1
  )

  # check that rownames of data.frame are used as default
  g <- glm(y ~ a, data=d)
  testthat::expect_identical(
              rownames(d),
              index(estimate(g))
            )

  # check that order is preserved (not sorted by default)
  d2 <- d[n:1,]
  g2 <- glm(y ~ a, data=d2)
  testthat::expect_identical(
              rev(rownames(d)),
              index(estimate(g2))
            )

  # user supplied id
  testthat::expect_identical(
              d$id,
              index(estimate(g, id=d$id))
            )

  # user supplied id, order preserved
  testthat::expect_identical(
              d$id2,
              index(estimate(g, id=d$id2))
            )

  # check that sort argument works
  testthat::expect_identical(
              sort(rownames(d)),
              index(merge(estimate(g2), sort=TRUE))
            )

  # check that id also works with transformations
  e <- estimate(g,
           function(x) x,
           id = d$id2)
  testthat::expect_identical(
              d$id2,
              index(e)
            )
  e1 <- estimate(g, id=d$id)
  testthat::expect_equivalent(
              IC(e),
              IC(e1),
              )
})
