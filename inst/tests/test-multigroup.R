context("Multiple Group")

test_that("Multiple group I", {
  m <- lvm(y~x)
  set.seed(1)
  d <- sim(m,100)
  ## Just a stratified analysis
  e <- estimate(list("Group A"=m,"Group B"=m),list(d,d))
  expect_equivalent(coef(e)[c(1,3)],coef(lm(y~x,d)))
  expect_equivalent(coef(e)[c(2,5)],coef(lm(y~x,d)))
})

test_that("Multiple group II", {
  m <- baptize(lvm(y~x))
  set.seed(1)
  d <- sim(m,100)
  ## Just a standard linear regression (single group)
  e <- estimate(list(m,m),list(d,d))
  expect_identical(coef(e,level=2)[[1]],coef(e,level=2)[[2]])
  expect_equivalent(coef(e,level=2)[[1]][1:2,1],coef(lm(y~x,cbind(d,d)))) 
})


context("Missing data")

test_that("Missing data analysis", {
  ## Random intercept model
  m <- lvm(c(y1,y2,y3)~x+u); latent(m) <- ~u
  set.seed(1)
  ## Missing on first two outcomes
  d <- makemissing(sim(m,200),p=0.3,cols=c("y1","y2"))  
  e <- estimate(m,d,missing=TRUE)
  expect_true("lvm.missing"%in%class(e))
  expect_true(sum(unlist(lapply(e$estimate$model$data,nrow)))==200)
  ## Convergence:
  g <- gof(e)
  expect_true(mean(score(e))<1e-3)
  expect_true(g$rankV==length(pars(e)))
})

test_that("Multiple group, missing data analysis", {
  m <- lvm(list(c(y1,y2,y3)~u,u~x)); latent(m) <- ~u
  m <- baptize(fixsome(m))
  regression(m,u~x) <- NA
  covariance(m,~u) <- NA
  set.seed(1)
  ## Missing on all outcomes
  d1 <- makemissing(sim(m,500),cols=c("y1","y2"),p=0.3)
  d2 <- makemissing(sim(m,500),cols=c("y1","y2"),p=0.3)
  e <- estimate(list(m,m),list(d1,d2),missing=TRUE)
  g <- gof(e)
  expect_true(g$n==1000)
  expect_true(mean(score(e))<1e-3)
  expect_true(g$rankV==length(pars(e)))
})






