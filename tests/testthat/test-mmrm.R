skip_if_not_installed("mmrm")
skip_if_not_installed("mets")
library("mmrm")

data(fev_data, package = "mmrm")

fevd <- subset(fev_data, select = -c(VISITN2, AVISIT)) |>
  mets::drename(FEV ~ FEV1)
dw <- mets::fast.reshape(fevd, id = "USUBJID",
                         varying = "FEV", num = "VISITN") |>
  subset(select = -c(FEV4))
dw <- na.omit(dw) # Data with complete-case in wide-format

dl <- mets::fast.reshape(
  dw,
  varying = list(FEV = c("FEV1", "FEV2", "FEV3")),
) |>
  transform(AVISIT = factor(num),
            ARMCD = factor(ARMCD),
            USUBJID = factor(USUBJID)) # Data in long format

ctrl <- list(eval.max = 1000, iter.max = 1000, rel.tol = 1e-9)
fit_mmrm <- function(formula, data = dl, ...) {
  mmrm(formula, data = data, reml = FALSE,
       optimizer = "nlminb", optimizer_control = ctrl, ...)
}

fits <- list(
  ## Unstructured covariance
  us = fit_mmrm(FEV ~ -1 + AVISIT + ARMCD:AVISIT + us(AVISIT | USUBJID)),
  ## Grouped unstructured covariance
  us_grp = fit_mmrm(FEV ~ -1 + AVISIT + ARMCD:AVISIT +
                      us(AVISIT | ARMCD / USUBJID)),
  ## Ante-dependence heterogeneous covariance
  adh = fit_mmrm(FEV ~ -1 + AVISIT + ARMCD:AVISIT + adh(AVISIT | USUBJID)),
  ## Ante-dependence homogeneous covariance
  ad = fit_mmrm(FEV ~ -1 + AVISIT + ARMCD:AVISIT + ad(AVISIT | USUBJID)),
  ## Auto-regressive covariance
  ar1h = fit_mmrm(FEV ~ -1 + AVISIT + ARMCD:AVISIT + ar1h(AVISIT | USUBJID)),
  ## Weighted model, unstructured covariance
  us_w = fit_mmrm(FEV1 ~ ARMCD * AVISIT + us(AVISIT | USUBJID),
                  data = fev_data, weights = fev_data$WEIGHT),
  ## Weighted model, toeplitz covariance structure
  toep_grp_w = fit_mmrm(FEV1 ~ ARMCD * AVISIT + toep(AVISIT | ARMCD / USUBJID),
                        data = fev_data, weights = fev_data$WEIGHT)
)

## log-likelihood as a function of parameters (mean (beta) and covariance
## (theta)).
.mmrm_loglik <- function(fit, beta = NULL, theta = NULL) {
  subj <- .mmrm_subjects(fit, theta = theta)
  if (is.null(beta)) beta <- mmrm::component(fit, "beta_est")
  ## Weighted covariance of a subject. mmrm models the covariance of subject i as
  ##   Sigma_i = W_i^{-1/2} Sigma W_i^{-1/2}, W_i = diag(w_i) i.e., cov(y_ij,
  ##   y_ik) = Sigma_jk / sqrt(w_ij * w_ik).
  Sinv <- lapply(subj, function(s) Inverse(.mmrm_wsigma(s$Sigma, s$w)))
  Sdet <- lapply(Sinv, function(s) attributes(s)$det)
  res <- lapply(subj, function(s) s$y - as.numeric(s$X %*% beta))
  loglik <- unlist(Map(function(Si, D, r) {
    -0.5 * (
      ncol(Si) * log(2 * pi) +
        log(D) +
        as.numeric(t(r) %*% Si %*% r)
    )
  }, Sinv, Sdet, res))
  sum(loglik)
}

test_that(".mmrm_subjects ordering of subjects and visits", {
  ref <- .mmrm_subjects(fit = fits$us, theta = NULL)

  new_dl <- dl
  new_dl$AVISIT <- factor(new_dl$AVISIT, levels = c("FEV3", "FEV1", "FEV2"))
  new_dl$USUBJID <- factor(
    new_dl$USUBJID,
    levels = levels(new_dl$USUBJID)[c(2, 1, 3:nlevels(new_dl$USUBJID))]
  )
  new_fit <- fit_mmrm(
    FEV ~ -1 + AVISIT + ARMCD:AVISIT + us(AVISIT | USUBJID),
    data = new_dl
  )
  ms <- .mmrm_subjects(fit = new_fit, theta = NULL)

  ## full frame in the model fit is ordered by id and visit levels
  new_ord <- c(2, 3, 1)
  expect_equal(ref[[1]]$visits, ms[[2]]$visits[new_ord])
  expect_equal(ref[[1]]$Sigma, ms[[2]]$Sigma[new_ord, new_ord],
               tolerance = 1e-6, ignore_attr = TRUE)
})

## test variance covariance constructor and parameter shifting (for numerical
## derivation)
test_that(".mmrm_varcor", {
  for (nam in names(fits)) {
    fit <- fits[[nam]]
    theta <- mmrm::component(fit, "theta_est")
    ref <- mmrm::VarCorr(fit)
    expect_equal(.mmrm_varcor(fit = fit, theta = theta), ref, label = nam)
    tmp <- .mmrm_varcor(fit = fit, theta = theta + 0.1)
    expect_true(all(abs(unlist(ref) - unlist(tmp)) > 1e-2), label = nam)
  }
})

## Weights (cov(y_ij, y_ik) = Sigma_jk / sqrt(w_ij * w_ik))
test_that(".mmrm_wsigma", {
  set.seed(1)
  S <- crossprod(matrix(rnorm(16), 4))
  w <- runif(4, 0.5, 2)
  D <- diag(1 / sqrt(w))
  expect_equal(.mmrm_wsigma(S, w), D %*% S %*% D, ignore_attr = TRUE)
  expect_identical(.mmrm_wsigma(S, rep(1, 4)), S)
})

test_that("score.mmrm", {
  for (nam in names(fits)) {
    fit <- fits[[nam]]
    ## score is zero at the MLE
    expect_lt(max(abs(colMeans(score(fit)))), 1e-6, label = nam)
    expect_equal(.mmrm_loglik(fit, beta = coef(fit)), as.numeric(logLik(fit)),
                 label = nam)
    ## Compare with numerical derivatives away from the MLE
    S0 <- numDeriv::jacobian(\(p) .mmrm_loglik(fit, beta = p), coef(fit) + 1)
    S <- .mmrm_score_beta(fit, beta = coef(fit) + 1)
    expect_equal(colSums(S), as.numeric(S0), tolerance = 1e-6,
                 ignore_attr = TRUE, label = nam)
    U0 <- numDeriv::jacobian(\(p) .mmrm_loglik(fit, theta = p),
                             fit$theta_est + 1)
    U <- .mmrm_score_theta(fit, theta = fit$theta_est + 1)
    expect_equal(colSums(U), as.numeric(U0), tolerance = 1e-6,
                 ignore_attr = TRUE, label = nam)
  }
})

## Transformation to real variance/covariance scale
test_that("mmrm2sigma, comparison with lava (unstructured covariance)", {
  fit <- fits$us
  f <- mmrm2sigma(fit)
  ea <- estimate(fit)
  avar <- transform(ea, f)
  # check id is there
  expect_equal(index(avar), sort(levels(dl$USUBJID)), ignore_attr = TRUE)
  expect_equal(dim(f(vec = FALSE)), c(3L, 3L))

  m <- lvm(c(FEV1, FEV2, FEV3) ~ ARMCD) |>
    covariance(~ FEV1 + FEV2 + FEV3, pairwise = TRUE)
  el <- estimate(estimate(m, dw))

  meanpar_mmrm <- subset(ea, 1:6)
  meanpar_lava <- subset(el, 1:6)
  varpar_mmrm <- subset(avar)
  varpar_lava <- subset(el, c(7, 10, 8, 11, 12, 9))

  expect_lt(mean((vcov(meanpar_mmrm) - vcov(meanpar_lava))^2), 1e-3)
  expect_lt(mean((coef(meanpar_mmrm) - coef(meanpar_lava))^2), 1e-9)
  expect_lt(mean((vcov(varpar_mmrm) - vcov(varpar_lava))^2), 1e-3)
  expect_lt(mean((coef(varpar_mmrm) - coef(varpar_lava))^2), 1e-9)
})

## visitn indexes the observed visits in the full covariance matrix
check_mmrm_visitn <- function(fit) {
  subj <- .mmrm_subjects(fit)
  Sig <- mmrm::VarCorr(fit)
  grouped <- is.list(Sig) && !is.matrix(Sig)
  for (s in subj) {
    Sfull <- if (grouped) Sig[[s$group]] else Sig
    if (grouped) {
      expect_true(s$group %in% names(Sig))
    } else {
      expect_true(is.na(s$group))
    }
    expect_type(s$visitn, "integer")
    expect_equal(s$visitn, match(s$visits, rownames(Sfull)))
    expect_equal(rownames(Sfull)[s$visitn], s$visits)
    expect_equal(Sfull[s$visitn, s$visitn, drop = FALSE], s$Sigma,
                 ignore_attr = TRUE)
  }
  ## Data contains subjects with missing visits
  expect_true(any(vapply(subj, function(s) length(s$visitn), 1L) <
                  ncol(if (grouped) Sig[[1]] else Sig)))
}

test_that(".mmrm_subjects with missing data", {
  ## Grouped toeplitz
  check_mmrm_visitn(fit_mmrm(
    FEV1 ~ -1 + AVISIT + ARMCD:AVISIT + toep(AVISIT | ARMCD / USUBJID),
    data = fev_data
  ))
  ## Toeplitz
  check_mmrm_visitn(fit_mmrm(
    FEV1 ~ -1 + AVISIT + ARMCD:AVISIT + toep(AVISIT | USUBJID),
    data = fev_data
  ))
})
