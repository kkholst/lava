test_that("score.glm can handle a overparameterized model (NA coef)", {
  set.seed(1)
  n <- 100
  data <- data.frame(
    id = seq_len(n),
    a = rbinom(n, 1, 0.5),
    x = rnorm(n)
  )
  data$x_dup <- data$x  # exact collinearity with `x`
  data$y <- 1 + data$a + data$x + rnorm(n)

  m_ref <- glm(y ~ a * x, data = data)

  score_ref <- score.glm(m_ref, indiv = FALSE)

  m <- glm(y ~ a * x + x_dup, data = data)
  score_na <- expect_warning(
    score.glm(m, indiv = FALSE),
    "Over-parameterized"
  )

  expect_equal(score_ref, score_na[!is.na(coef(m))])

  expect_true(score_na["x_dup"] == 0)

  score_ref <- score.glm(m_ref, indiv = TRUE)
  score_na <- suppressWarnings(score.glm(m, indiv = TRUE))

  expect_equal(score_ref, score_na[,!is.na(coef(m))],
               check.attributes = FALSE)

  expect_equal(
    attr(score_ref, "bread"),
    attr(score_na, "bread")[names(coef(m)) != "x_dup",
                            names(coef(m)) != "x_dup"]
  )

  expect_true(
    all(attr(score_na, "bread")[names(coef(m)) == "x_dup", ] == 0) &&
      all(attr(score_na, "bread")[, names(coef(m)) == "x_dup"] == 0)
  )

  ## check IC

  expect_equal(
    IC(m_ref),
    suppressWarnings(IC(m))[, names(coef(m)) != "x_dup"],
    check.attributes = FALSE
  )

})
