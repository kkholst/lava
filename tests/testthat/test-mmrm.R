
library("mmrm")
library("lava")

data(fev_data, package="mmrm")

fevd <- subset(fev_data, select=-c(VISITN2, AVISIT)) |>
  mets::drename(FEV ~ FEV1)
dw <- mets::fast.reshape(fevd, id="USUBJID",
                         varying="FEV", num="VISITN") |>
  subset(select=-c(FEV4))
dw <- na.omit(dw) # Data with complete-case in wide-format

dl <- mets::fast.reshape(
              dw,
              varying = list(FEV=c("FEV1", "FEV2", "FEV3")),
              ) |>
  transform(AVISIT=factor(num),
            ARMCD=factor(ARMCD),
            USUBJID=factor(USUBJID)) # Data in long format

ctrl <- list(eval.max = 1000, iter.max = 1000, rel.tol = 1e-9)
## Unstructured covariance
a <- mmrm(
  FEV ~ -1 + AVISIT + ARMCD : AVISIT + us(AVISIT | USUBJID),
  data = dl, reml = FALSE,
  optimizer = "nlminb",
  optimizer_control = ctrl
)

## Grouped unstructured covariance
b <- mmrm(FEV ~ -1 + AVISIT + ARMCD:AVISIT + us(AVISIT | ARMCD / USUBJID),
          data = dl,
          reml = FALSE,
          optimizer = "nlminb",
          optimizer_control = ctrl)

## Ante-dependence heterogeneous covariance
c <- mmrm(
  FEV ~ -1 + AVISIT + ARMCD : AVISIT + adh(AVISIT | USUBJID),
  data = dl, reml = FALSE,
  optimizer = "nlminb",
  optimizer_control = ctrl
)

## Ante-dependence homogeneous covariance
d <- mmrm(
  FEV ~ -1 + AVISIT + ARMCD : AVISIT + ad(AVISIT | USUBJID),
  data = dl, reml = FALSE,
  optimizer = "nlminb",
  optimizer_control = ctrl
)

## Auto-regressive covariance
e <-  mmrm(
  FEV ~ -1 + AVISIT + ARMCD : AVISIT + ar1h(AVISIT | USUBJID),
  data = dl, reml = FALSE,
  optimizer = "nlminb",
  optimizer_control = ctrl
)

## Weighted model, unstructured covariance
aw <- mmrm(FEV1 ~ ARMCD * AVISIT + us(AVISIT | USUBJID),
           data = fev_data, reml = FALSE, weights = fev_data$WEIGHT,
           optimizer = "nlminb", optimizer_control = ctrl)

## Weighted model, toeplitz covariance structure
bw <- mmrm(FEV1 ~ ARMCD * AVISIT + toep(AVISIT | ARMCD / USUBJID),
           data = fev_data, reml = FALSE, weights = fev_data$WEIGHT,
           optimizer = "nlminb", optimizer_control = ctrl)

## log-likelihood as a function of parameters (mean (beta) and covariance
## (theta)).
.mmrm_loglik <- function(fit, beta=NULL, theta=NULL) {
  subj <- targeted:::.mmrm_subjects(fit, theta=theta)
  if (is.null(beta)) beta <- mmrm::component(fit, "beta_est")
  ## Weighted covariance of a subject. mmrm models the covariance of subject i as
  ##   Sigma_i = W_i^{-1/2} Sigma W_i^{-1/2}, W_i = diag(w_i) i.e., cov(y_ij,
  ##   y_ik) = Sigma_jk / sqrt(w_ij * w_ik).
  Sinv <- lapply(subj, function(s) lava::Inverse(targeted:::.mmrm_wsigma(s$Sigma, s$w)))
  Sdet <- lapply(Sinv, function(s) attributes(s)$det)
  res    <- lapply(subj, function(s) s$y - as.numeric(s$X %*% beta))
  loglik <- unlist(Map(function(Si, D, r) {
    -0.5 * (
      ncol(Si) * log(2 * pi) +
        log(D) +
        as.numeric(t(r) %*% Si %*% r)
    )
  }, Sinv, Sdet, res))
  sum(loglik)
}

## test ordering of subjects in .mmrm_subjects
test_mmrm_subject <- function(fit, data_long) {
  ref <- targeted:::.mmrm_subjects(fit = fit, theta = NULL)

  new_data_long <- data_long
  new_data_long$AVISIT <- factor(new_data_long$AVISIT, levels = c("FEV3", "FEV1", "FEV2"))
  new_data_long$USUBJID <- factor(new_data_long$USUBJID,
                                  levels = levels(new_data_long$USUBJID)[c(2, 1, 3:65)])

  new_fit <- mmrm(
    FEV ~ -1 + AVISIT + ARMCD : AVISIT + us(AVISIT | USUBJID),
    data = new_data_long, reml = FALSE,
    optimizer = "nlminb",
    optimizer_control = list(
      eval.max = 1000,
      iter.max = 1000,
      rel.tol = 1e-9
    )
  )

  ms <- targeted:::.mmrm_subjects(fit = new_fit, theta = NULL)

  ## full frame in the model fit is ordered by id and visit levels
  new_ord <- c(2, 3, 1)
  expect_true(
    all(ref[[1]]$visits == ms[[2]]$visits[new_ord])
  )
  expect_true(
    all(ref[[1]]$Sigma - ms[[2]]$Sigma[new_ord, new_ord] < 10-8)
  )

}

test_mmrm_subject(fit = a, data_long = dl)

## test variance covariance constructor
## test parameter shifting (for numerical derivation)
test_mmrm_varcor <- function(fit) {
  ref <- mmrm::VarCorr(fit)
  tmp <- targeted:::.mmrm_varcor(
                      fit = fit,
                      theta = mmrm::component(fit, "theta_est")
                    )
  expect_equal(
    ref,
    tmp
  )

  tmp <- targeted:::.mmrm_varcor(
                      fit = fit,
                      theta = mmrm::component(fit, "theta_est") + 0.1
                    )
  if (is.list(ref)) {
    expect_true(
      all(abs(unlist(ref) - unlist(tmp)) > 1e-2)
    )
  } else {
    expect_true(
      all(abs(ref - tmp) > 1e-2)
    )
  }


}
test_mmrm_varcor(a)
test_mmrm_varcor(b)
test_mmrm_varcor(c)
test_mmrm_varcor(d)
test_mmrm_varcor(e)
test_mmrm_varcor(aw)
test_mmrm_varcor(bw)

## Weights (cov(y_ij, y_ik) = Sigma_jk / sqrt(w_ij * w_ik))
test_mmrm_wsigma <- function() {
  set.seed(1)
  S <- crossprod(matrix(rnorm(16), 4))
  w <- runif(4, 0.5, 2)
  D <- diag(1 / sqrt(w))
  expect_equivalent(targeted:::.mmrm_wsigma(S, w), D %*% S %*% D)
  expect_identical(targeted:::.mmrm_wsigma(S, rep(1, 4)), S)
}
test_mmrm_wsigma()

## test score function
test_mmrm_score <- function(fit) {
  opt <- lava.options(Dmethod="Richardson")
  expect_true(mean(colMeans(score(fit)))<1e-9)
  ll <- c(.mmrm_loglik(fit, beta=coef(fit)), logLik(fit))
  S0 <- numDeriv::jacobian(\(p) .mmrm_loglik(fit, beta=p),
                           coef(fit)+1)
  S <- targeted:::.mmrm_score_beta(fit, beta=coef(fit)+1)
  U0 <- numDeriv::jacobian(\(p) .mmrm_loglik(fit, theta=p), fit$theta_est+1)
  U <- targeted:::.mmrm_score_theta(fit, theta=fit$theta_est+1)

  expect_equivalent(ll[1], ll[2])
  expect_true(mean(as.numeric(S0)-colSums(S))<1e-6)
  expect_true(mean(as.numeric(U0)-colSums(U))<1e-6)
  lava.options(opt)
}
test_mmrm_score(a)
test_mmrm_score(b)
test_mmrm_score(c)
test_mmrm_score(d)
test_mmrm_score(e)
test_mmrm_score(aw)
test_mmrm_score(bw)

# Transformation to real variance/covariance scale
## lava benchmark, unstructured covariance
test_mmrm_sigma <- function() {
  f <- targeted:::mmrm2sigma(a)
  ea <- estimate(a)
  avar <- transform(ea, f)
  # check id is there
  expect_equivalent(
              index(avar), sort(levels(dl$USUBJID))
  )
  covarest <- f(vec=FALSE)
  expect_true(all(dim(covarest) == c(3L, 3L)))

  ## Comparison with lava
  m <- lvm(c(FEV1,FEV2,FEV3) ~ ARMCD) |>
    covariance(~FEV1+FEV2+FEV3, pairwise=TRUE)
  es <- estimate(m, dw)
  score(es)
  e <- estimate(es)

  meanpar_mmrm <- subset(ea, 1:6)
  meanpar_lava <- subset(e, 1:6)
  varpar_mmrm <- subset(avar)
  varpar_lava <- subset(e, c(7,10,8,11,12,9))

  expect_true(
    mean((vcov(meanpar_mmrm)-vcov(meanpar_lava))^2)<1e-3
  )
  expect_true(
    mean((coef(meanpar_mmrm)-coef(meanpar_lava))^2)<1e-9
  )
  expect_true(
    mean((vcov(varpar_mmrm)-vcov(varpar_lava))^2)<1e-3
  )
  expect_true(
    mean((coef(varpar_mmrm)-coef(varpar_lava))^2)<1e-9
  )
}
test_mmrm_sigma()

## Missing data, grouped toeplitz
a <- mmrm(
  FEV1 ~ -1 + AVISIT + ARMCD : AVISIT + toep(AVISIT | ARMCD / USUBJID),
  data = fev_data, reml = FALSE,
  optimizer = "nlminb",
  optimizer_control = list(
    eval.max = 1000,
    iter.max = 1000,
    rel.tol = 1e-9
  )
)

## visitn indexes the observed visits in the full covariance matrix
test_mmrm_visitn <- function(fit) {
  subj <- targeted:::.mmrm_subjects(fit)
  Sig <- mmrm::VarCorr(fit)
  grouped <- is.list(Sig) && !is.matrix(Sig)
  for (s in subj) {
    Sfull <- if (grouped) Sig[[s$group]] else Sig
    if (grouped) {
      expect_true(s$group %in% names(Sig))
    } else {
      expect_true(is.na(s$group))
    }
    expect_true(is.integer(s$visitn))
    expect_equal(s$visitn, match(s$visits, rownames(Sfull)))
    expect_equal(rownames(Sfull)[s$visitn], s$visits)
    expect_equivalent(Sfull[s$visitn, s$visitn, drop = FALSE], s$Sigma)
  }
  ## Data contains subjects with missing visits
  expect_true(any(vapply(subj, function(s) length(s$visitn), 1L) <
                  ncol(if (grouped) Sig[[1]] else Sig)))
}
test_mmrm_visitn(a)

## Missing data, toeplitz
a <- mmrm(
  FEV1 ~ -1 + AVISIT + ARMCD : AVISIT + toep(AVISIT | USUBJID),
  data = fev_data, reml = FALSE,
  optimizer = "nlminb",
  optimizer_control = list(
    eval.max = 1000,
    iter.max = 1000,
    rel.tol = 1e-9
  )
)
test_mmrm_visitn(a)
