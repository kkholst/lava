library("mmrm")
library("lava")

data(fev_data, package = "mmrm")
ctrl <- list(eval.max = 1000, iter.max = 1000, rel.tol = 1e-9)

## Lightweight estimate.mmrm object on the sigma scale (avoids computing the
## influence function)
est_obj <- function(fit) {
  structure(list(coef = c(fit$beta_est, targeted:::mmrm2sigma(fit)()), fit = fit,
                 sigma = TRUE),
            class = c("estimate.mmrm", "estimate"))
}

## mmrm's conditional prediction at visit 'time' for each subject
pred_mmrm <- function(fit, data, time) {
  fp <- fit$formula_parts
  d <- complete_visits(data, fp$subject_var, fp$visit_var, fp$response_var,
                   levels(fit$tmb_data$full_frame[[fp$visit_var]]))
  pm <- predict(fit, newdata = d, conditional = TRUE)
  rr <- which(as.integer(d[[fp$visit_var]]) == time)
  setNames(pm[rr], as.character(d[[fp$subject_var]][rr]))
}

## complete_visits
test_complete_visits <- function() {
  d <- data.frame(id = c("b", "b", "a", "a", "c"),
                  visit = c(3, 1, 1, 2, 2),
                  y = c(1, 2, 3, 4, 5), x = c(10, 20, 30, 40, 50))
  res <- complete_visits(d, "id", "visit", "y")
  expect_equal(nrow(res), 9L)
  expect_equal(res$id, rep(c("b", "a", "c"), each = 3))
  expect_equal(res$visit, factor(rep(1:3, 3)))
  expect_equal(res$y, c(2, NA, 1, 3, 4, NA, NA, 5, NA))
  ## Added rows: previous row (LOCF), or next row if first visit missing
  expect_equal(res$x, c(20, 20, 10, 30, 40, 40, 50, 50, 50))
  ## Without response: added rows keep the copied response
  expect_equal(complete_visits(d, "id", "visit")$y[2], 2)
  ## Visit levels
  res <- complete_visits(d, "id", "visit", "y", levels = 1:4)
  expect_equal(nlevels(res$visit), 4L)
  expect_equal(nrow(res), 12L)
  ## Complete data is unchanged (up to ordering)
  res2 <- complete_visits(res, "id", "visit", "y")
  expect_equal(res2, res)
  ## Errors
  expect_error(complete_visits(d, "id", "visit", levels = 1:2))
  expect_error(complete_visits(rbind(d, d[1, ]), "id", "visit"))
}
test_complete_visits()

## Comparison with the conditional predictions of mmrm
test_predict_vs_mmrm <- function(formula) {
  ## Data with entire rows removed
  dsub <- subset(fev_data, !is.na(FEV1) | VISITN %% 2 == 0)
  for (w in list(NULL, fev_data$WEIGHT)) {
    fit <- mmrm(formula, data = fev_data, reml = FALSE, weights = w,
                optimizer = "nlminb", optimizer_control = ctrl)
    e <- est_obj(fit)
    for (time in 1:4) {
      pr <- predict(e, time = time)
      expect_equal(length(pr), nlevels(fev_data$USUBJID))
      expect_equivalent(pr, pred_mmrm(fit, fev_data, time)[names(pr)],
                        tolerance = 1e-8)
      obs <- attr(pr, "observed")
      y <- fev_data$FEV1[fev_data$AVISIT == levels(fev_data$AVISIT)[time]]
      expect_equal(obs, !is.na(y))
      expect_equivalent(pr[obs], y[obs])

      pr <- predict(e, newdata = dsub, time = time)
      expect_equivalent(pr, pred_mmrm(fit, dsub, time)[names(pr)],
                        tolerance = 1e-8)
    }
  }
}
test_predict_vs_mmrm(FEV1 ~ ARMCD * AVISIT + us(AVISIT | USUBJID))
test_predict_vs_mmrm(FEV1 ~ ARMCD * AVISIT + toep(AVISIT | ARMCD / USUBJID))


## Delta method
test_predict_delta <- function() {
  d <- subset(fev_data, USUBJID %in% levels(USUBJID)[1:80])
  fit <- mmrm(FEV1 ~ ARMCD * AVISIT + us(AVISIT | USUBJID), data = d,
              reml = FALSE, optimizer = "nlminb", optimizer_control = ctrl)
  e <- estimate(fit, sigma = TRUE)
  newd <- subset(d, USUBJID %in% c("PT1", "PT2", "PT3"))
  ## PT1 is unobserved at the first visit
  pr <- predict(e, newdata = newd, time = 1)
  expect_true(!all(attr(pr, "observed")))
  ep <- estimate(e, function(p) {
    mean(predict(e, p = p, newdata = newd, time = 1))
  })
  expect_equivalent(coef(ep), mean(pr))
  expect_true(is.finite(vcov(ep)) && vcov(ep) > 0)
  ## Only observed outcomes: no uncertainty
  expect_true(all(attr(predict(e, newdata = newd), "observed")))
  ep <- estimate(e, function(p) mean(predict(e, p = p, newdata = newd)))
  expect_equivalent(as.numeric(vcov(ep)), 0)
}
test_predict_delta()
