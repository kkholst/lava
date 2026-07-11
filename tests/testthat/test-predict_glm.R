context("predict_glm")

set.seed(42)
n <- 100
d <- data.frame(y = rpois(100, lambda = 3), x1 = rnorm(100), x2 = rnorm(100))

test_that("predict_glm matches predict.glm for poisson(log)", {
  m <- glm(y ~ x1 + x2, data = d, family = poisson(link = "log"))
  pred_resp <- as.numeric(predict_glm(m, type = "response"))
  pred_link <- as.numeric(predict_glm(m, type = "link"))
  expect_equal(pred_resp, as.numeric(predict(m, type = "response")))
  expect_equal(pred_link, as.numeric(predict(m, type = "link")))

  ## with modified parameter
  p0 <- coef(m) * 0
  pred_mod <- as.numeric(predict_glm(m, p = p0, type = "response"))
  expect_equal(pred_mod, exp(rep(0, nrow(d))))

  # with data
  pr <- as.numeric(predict_glm(m, p = p0, data = head(d, 5L)))
  expect_true(length(pr) == 5L)
  expect_true(all(pr == 1L))

  # with user offset
  pr <- as.numeric(predict_glm(m, p = p0, data = head(d, 5L), offset = 1L))
  expect_true(all(pr == exp(1)))
})

test_that("predict_glm matches predict.glm with offset (poisson)", {
  m <- glm(y ~ x1 + offset(x2), data = d, family = poisson())
  pred_resp <- as.numeric(predict_glm(m, type = "response"))
  pred_link <- as.numeric(predict_glm(m, type = "link"))
  expect_equal(pred_resp, as.numeric(predict(m, type = "response")))
  expect_equal(pred_link, as.numeric(predict(m, type = "link")))
})
