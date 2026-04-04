context("composite likelihood")


test_that("complik", {
  m <- lvm(c(y1,y2,y3)~b*x+1*u[0],latent=~u)
  ordinal(m,K=2) <- ~y1+y2+y3
  d <- sim(m,50,seed=1)
  if (requireNamespace("mets", quietly=TRUE)) {
    e1 <- complik(m,d,control=list(trace=0),type="all")
    dl <- mets::fast.reshape(d)
    g <- lme4::glmer(y~x+(1|id), family=binomial("probit"), data=dl)
    # MLE and complik should agree reasonably
    # even though e1 is fitted with different intercepts so the
    # two models are not identical (lmer nested within complik)
    expect_true((coef(e1)["y1~x"]-lme4::fixef(g)["x"])<0.1)
    expect_true((vcov(e1)["y1~x","y1~x"]-vcov(g)["x","x"])<0.025)

    expect_equivalent(vcov(e1), var_ic(IC(e1)))
    expect_true(inherits(logLik(e1), "logLik"))
    expect_true(mean(score(e1)^2)<1e-4)
  }
})
