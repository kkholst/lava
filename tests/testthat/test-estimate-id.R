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

test_that("estimate index from rownames and `id` arg", {
  ic0 <- ic1 <- cbind(1:5-3)
  rownames(ic1) <- paste0("i", 1:5)
  e0 <- estimate(coef=0, IC=ic0) # no index
  expect_true(is.null(index(e0)))
  e1 <- estimate(coef=0, IC=ic1)
  expect_identical(rownames(ic1), index(e1))
  e <- estimate(e0, id=1:5)
  expect_identical(rownames(IC(e)), as.character(1:5))
})

test_that("id as column name of data", {
  set.seed(1)
  d <- data.frame(x = rnorm(40))
  d$y <- d$x + rnorm(40)
  d$cl <- rep(paste0("c", 1:10), each = 4)
  g <- lm(y ~ x, data = d)
  # clustering (no averaging)
  e <- estimate(g, id = "cl", data = d)
  expect_equal(vcov(e), vcov(estimate(g, id = d$cl)))
  expect_identical(index(e), unique(d$cl))
  # standardization: id refers to rows of 'data'
  f <- function(p, data) p[1] + p[2] * data[, "x"]
  ec <- estimate(g, id = d$cl)
  a <- estimate(ec, f, data = d, id = "cl", average = TRUE)
  expect_equal(
    vcov(a),
    vcov(estimate(ec, f, data = d, id = d$cl, average = TRUE))
  )
  expect_identical(index(a), unique(d$cl))
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

test_that("id=NULL removes the id (index)", {
  set.seed(1)
  n <- 40
  d <- data.frame(y = rnorm(n), x = rnorm(n))
  rownames(d) <- paste0("r", seq_len(n))
  d$id <- paste0("a", seq_len(n))
  d$cl <- rep(paste0("c", 1:10), each = 4)
  g <- glm(y ~ x, data = d)
  e <- estimate(g)
  e0 <- estimate(g, id = NULL)
  expect_null(index(e0))
  expect_null(rownames(IC(e0)))
  expect_equal(vcov(e0), vcov(e))
  expect_equal(coef(e0), coef(e))
  expect_equivalent(IC(e0), IC(e))

  # estimate objects and index<-
  e1 <- estimate(estimate(g, id = d$id), id = NULL)
  expect_null(index(e1))
  expect_null(rownames(IC(e1)))
  expect_equal(vcov(e1), vcov(e))
  e2 <- e
  index(e2) <- NULL
  expect_null(index(e2))
  expect_null(rownames(IC(e2)))

  # clustered IF is kept, only the ids are removed
  ec <- estimate(g, id = d$cl)
  ec0 <- estimate(ec, id = NULL)
  expect_null(index(ec0))
  expect_equal(vcov(ec0), vcov(ec))
  expect_equivalent(IC(ec0), IC(ec))

  # averaging: default linking, ids removed from result
  f <- function(p, data) p[1] + p[2] * data[, "x"]
  a <- estimate(g, f, data = d, average = TRUE)
  a0 <- estimate(g, f, data = d, id = NULL, average = TRUE)
  expect_null(index(a0))
  expect_null(rownames(IC(a0)))
  expect_equal(vcov(a0), vcov(a))

  # merge requires ids (or explicit independence / pairing)
  expect_error(merge(e0, e0), "Need id")
  expect_equal(vcov(merge(e0, e0, id = NULL))[1:2, 1:2], vcov(e))
  expect_equal(vcov(merge(e0, e0, paired = TRUE)), vcov(merge(e, e)))
})

test_that("cluster aggregation of IF (stack)", {
  set.seed(1)
  n <- 40
  d <- data.frame(y = rnorm(n), x = rnorm(n))
  d$cl <- rep(sample(paste0("c", 1:10)), each = 4) # unsorted cluster ids
  g <- glm(y ~ x, data = d)
  e <- estimate(g, id = d$cl)

  # first-appearance order and manual aggregation
  ic0 <- IC(g)
  ic1 <- rowsum(ic0, d$cl, reorder = FALSE) * 10 / n
  expect_identical(index(e), unique(d$cl))
  expect_equivalent(IC(e), ic1)
  expect_equal(rownames(IC(e)), unique(d$cl))

  # attributes are kept: 'bread' and number of observations 'N'
  expect_equivalent(attr(IC(e), "bread"), attr(ic0, "bread"))
  expect_equal(attr(IC(e), "N"), n)
  expect_no_error(summary(e, type = "mbn"))

  # re-clustering an estimate object keeps the original 'N'
  cl2 <- sub("c([0-9]+)", "k\\1", unique(d$cl))
  cl2[1:2] <- "k0"
  e2 <- estimate(e, id = cl2)
  expect_equal(attr(IC(e2), "N"), n)
  expect_equal(nrow(IC(e2)), 9)

  # missing values in id
  idna <- d$cl
  idna[3] <- NA
  expect_error(estimate(g, id = idna), "Missing values in 'id'")
})
