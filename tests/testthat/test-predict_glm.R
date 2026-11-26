# predict_glm

set.seed(42)
n <- 100
d <- data.frame(y = rpois(100, lambda = 3),
                x1 = rnorm(100), x2 = rnorm(100))

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

  # data with response variable
  pr <- as.numeric(predict_glm(m, p = p0, data=head(d, 5L)))
  expect_true(length(pr)==5L)
  expect_true(all(pr==1L))

  newd <- head(d, 5L)[, c("x1", "x2")]
  pr1 <- as.numeric(predict_glm(m, p = p0, data=newd))
  expect_equal(pr, pr1)

  # with user offset
  pr <- as.numeric(predict_glm(m, p = p0, data=newd, offset=1L))
  expect_true(all(pr==exp(1)))

  # formula environment is handled correctly in predict_glm when removing the
  # response variable
  m <- local({
    thr <- 2
    glm(y ~ x1 + I(x2 < thr), data = d, family = poisson(link = "log"))
  })
  pr <- as.numeric(predict_glm(m, data=newd))
  expect_equal(pr, as.numeric(predict(m, newd, type = "response")))
})

test_that("predict_glm matches predict.glm with offset (poisson)", {
  m <- glm(y ~ x1 + offset(x2), data = d, family = poisson())
  pred_resp <- as.numeric(predict_glm(m, type = "response"))
  pred_link <- as.numeric(predict_glm(m, type = "link"))
  expect_equal(pred_resp, as.numeric(predict(m, type = "response")))
  expect_equal(pred_link, as.numeric(predict(m, type = "link")))
})

test_that("predict_glm with offset and 'data' of different length", {
  d2 <- transform(d, off = x2)
  newd <- head(d2, 5L)

  ## offset() term of the formula
  m <- glm(y ~ x1 + offset(x2), data = d2, family = poisson())
  expect_equal(as.numeric(predict_glm(m, data = newd)),
               as.numeric(predict(m, newdata = newd, type = "response")))

  ## 'offset' argument of the glm() call
  m <- glm(y ~ x1, data = d2, family = poisson(), offset = off)
  expect_equal(as.numeric(predict_glm(m, data = newd)),
               as.numeric(predict(m, newdata = newd, type = "response")))

  ## A user supplied offset takes precedence
  m <- glm(y ~ x1, data = d2, family = poisson())
  expect_equal(as.numeric(predict_glm(m, p = coef(m) * 0, data = newd,
                                      offset = 1L)),
               rep(exp(1), 5L))
})

test_that("predict_glm handles factors with levels unused in the fit", {

  set.seed(1)
  n <- 200
  dd <- data.frame(
    t = factor(sample(c("A", "B", "C", "D"), n, TRUE),
                 levels = c("A", "B", "C", "D")),
    t2 = factor(sample(c("N", "Y"), n, TRUE), levels = c("N", "Y")),
    b = rnorm(n)
  )
  dd$y <- rnorm(n)

  ## Only two of the four t levels are present, but the levels attribute
  ## retains all four
  dds <- subset(dd, t %in% c("A", "B"))
  expect_length(levels(dds$t), 4L)

  m <- glm(y ~ t * t2 * b,
           data = dds)

  ## glm() drops the unused levels: fewer coefficients than a 4-level design
  expect_equal(m$xlevels$t, c("A", "B"))

  estobj <- estimate(m)
  pr <- predict(estobj, p = coef(estobj), data = dds)
  expect_equal(as.numeric(pr),
               as.numeric(predict(m, newdata = dds,
                                  type = "response")))
  expect_equal(NCOL(attr(pr, "grad")), length(coef(estobj)))

  ## Model fitted on a subset, predictions on data with all levels present
  expect_error(
    predict_glm(m, data = dd),
    "factor t has new levels"
  )
})

test_that("predict_glm respects factor level order and contrasts of the fit", {

  set.seed(1)
  n <- 100
  dd <- data.frame(x = rnorm(n),
                   g = factor(sample(c("a", "b", "c"), n, TRUE)))
  dd$y <- rnorm(n)
  newd <- head(dd, 5L)

  m <- glm(y ~ g + x, data = dd)
  ## Levels in a different order than in the fit
  newd2 <- transform(newd, g = factor(as.character(newd$g),
                                      levels = c("c", "b", "a")))
  expect_equal(as.numeric(predict_glm(m, data = newd2)),
               as.numeric(predict(m, newdata = newd2, type = "response")))

  ## Non-default contrasts
  m3 <- glm(y ~ g + x, data = dd, contrasts = list(g = "contr.sum"))
  expect_equal(as.numeric(predict_glm(m3, data = newd)),
               as.numeric(predict(m3, newdata = newd, type = "response")))

  ## A level that was not present in the fit is an error (as in predict.glm)
  newd3 <- transform(newd, g = factor(c("z", "a", "a", "b", "b")))
  expect_error(predict_glm(m3, data = newd3), "new level")

  ## Mismatching parameter vector gives an informative error
  expect_error(predict_glm(m, p = coef(m)[1:2], data = newd),
               "does not match the length")
})

test_that("predict_glm evaluates basis functions with the basis of the fit", {

  set.seed(1)
  n <- 100
  dd <- data.frame(x = rnorm(n))
  dd$y <- rnorm(n)
  newd <- head(dd, 5L)

  m <- glm(y ~ poly(x, 2), data = dd)
  expect_equal(as.numeric(predict_glm(m, data = newd)),
               as.numeric(predict(m, newdata = newd, type = "response")))
  ## Without 'data' the model matrix of the fit is used
  expect_equal(as.numeric(predict_glm(m)), as.numeric(fitted(m)))
  ## predict.estimate.glm must not materialize model.frame(fit) here
  e <- estimate(m)
  expect_equal(as.numeric(predict(e)), as.numeric(fitted(m)))
  expect_equal(as.numeric(predict(e, data = newd)),
               as.numeric(predict(m, newdata = newd, type = "response")))
})

test_that("score.glm/score.lm with 'data' use the design matrix of the fit", {

  set.seed(1)
  n <- 200
  dd <- data.frame(
    g = factor(sample(c("A", "B", "C", "D"), n, TRUE),
               levels = c("A", "B", "C", "D")),
    h = factor(sample(c("N", "Y"), n, TRUE), levels = c("N", "Y")),
    x = rnorm(n)
  )
  dd$y <- rnorm(n)
  dsub <- subset(dd, g %in% c("A", "B"))

  m <- glm(y ~ g * h * x, data = dsub)
  expect_equal(score.glm(m, data = dsub), score.glm(m))
  expect_true(max(abs(score.glm(m, data = dsub))) < 1e-8) # score at the MLE

  ml <- lm(y ~ g * h * x, data = dsub)
  expect_equal(score.lm(ml, data = dsub), score.lm(ml))

  ## Basis functions
  m2 <- glm(y ~ poly(x, 2), data = dsub)
  expect_equal(score.glm(m2, data = dsub), score.glm(m2))
})
