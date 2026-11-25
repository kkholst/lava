test_that("predict.estimate.glm accepts 'newdata' as alias for 'data'", {
  set.seed(1)
  d <- data.frame(x = rnorm(50), z = rnorm(50))
  d$y <- 1 + d$x + rnorm(50)
  m <- glm(y ~ x + z, data = d)
  e <- estimate(m)

  nd <- d[1:5, ]

  ## 'newdata' is honored (previously silently swallowed by '...')
  pr_newdata <- predict(e, newdata = nd)
  expect_equal(length(pr_newdata), nrow(nd))

  ## 'data' still works and matches 'newdata'
  pr_data <- predict(e, data = nd)
  expect_equal(pr_data, pr_newdata)

  ## Default (neither supplied) still equals original model.frame predictions
  pr_default <- predict(e)
  expect_equal(length(pr_default), nrow(d))
})

test_that("predict.estimate.glm errors on conflicting 'newdata' and 'data'", {
  set.seed(1)
  d <- data.frame(x = rnorm(30))
  d$y <- 1 + d$x + rnorm(30)
  m <- glm(y ~ x, data = d)
  e <- estimate(m)

  nd1 <- d[1:3, ]
  nd2 <- d[4:6, ]
  expect_error(predict(e, newdata = nd1, data = nd2), "not both")

  ## Identical values pass through silently
  expect_equal(
    predict(e, newdata = nd1, data = nd1),
    predict(e, newdata = nd1)
  )
})

test_that("expectation(predict, ...) still works for estimate.glm (vignette)", {
  set.seed(1)
  n <- 100
  d <- data.frame(id = seq_len(n), x = rnorm(n), a = rbinom(n, 1, 0.5))
  d$y <- plogis(-0.5 + d$a + d$x)
  d$y <- rbinom(n, 1, d$y)
  g <- glm(y ~ a + x, data = d, family = binomial)
  e <- estimate(g, id = d$id)

  ea <- expectation(predict, object = e, data = transform(d, a = 1),
                    id = "id", labels = "trt=1")
  expect_s3_class(ea, "estimate")
  expect_true(is.finite(coef(ea)))
})
