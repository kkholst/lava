context("ordreg: cumulative link regression")

test_that("ordreg",{
  set.seed(1)
  n <- 5e3
  x <- cbind(rnorm(n))
  J <- 5
  as <- c(seq(-1, 1, length.out=J-1))
  bs <- c(1)
  cumpr <- matrix(0, ncol=J-1, nrow=n)
  link <- pnorm # lava::expit
  for (j in seq_len(J-1)) {
    cumpr[, j] <- link(as[j] - bs*x)
  }
  # P(Y<=j|X=x) = expit(a_j + b*x)
  pr <- cbind(cumpr[,1], t(apply(cbind(cumpr, 1), 1, diff)))
  y <- apply(pr, 1, function(x) which(rmultinom(1, 1, x)==1L)-1)
  d <- data.frame(x=x, y=y)

  ## e <- ordinal::clm(ordered(y) ~ x, data=d)
  a <- ordreg(ordered(y) ~ x, data=d, family=binomial("probit"))
  ## estimate(a, function(x) lava:::ordreg_threshold(x[seq_len(J-1)]))
  s <- summary(a)
  cc <- coef(s$coef)
  expect_equivalent(cc, c(as, bs), tolerance=0.1)
}
)
