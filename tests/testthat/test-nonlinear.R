library(testthat)

context("nonlinear lvm")

## simulated data
f <- function(x) - 0.25*x^2
m <- lvm(x1+x2+x3~eta1, y1+y2+y3~eta2, latent=~eta1+eta2)
functional(m, eta2~eta1) <- f
d <- sim(m,1000,seed=1,latent=TRUE)

## Setup two-stage model
m1 <- lvm(x1+x2+x3~eta1,latent=~eta1)
m2 <- lvm(y1+y2+y3~eta2,latent=~eta2)
mm <- twostage(m1,m2,formula=eta2~eta1,type="spline")

test_that("twostage", {
  nonlinear(m2,type="quadratic") <- eta2~eta1
  a <- twostage(m1,m2,data=d)
  expect_true(abs(coef(a)["eta2~eta1_1"])<0.1)
  expect_true(abs(coef(a)["eta2~eta1_2"]+0.25)<0.1)
})

test_that("cv", {
  m <- list(lm(Sepal.Length~1, data=iris),
           lm(Sepal.Length~Species, data=iris),
           lm(Sepal.Length~Species * Petal.Length, data=iris))
  x <- cv(m, rep=2, nfolds=5, data=iris)
  expect_equivalent(dim(x$cv), c(2,5,3,1))
  rmse.1 <- mean(x$cv[,,1,1])
  rmse.2 <- mean(x$cv[,,2,1])
  rmse.3 <- mean(x$cv[,,3,1])
  expect_true(rmse.1>rmse.2)
  expect_true(rmse.2>rmse.3)
})

## testthat("twostage cv", {
##   ma <- lvm(c(x1,x2,x3)~u,latent=~u)
##   ms <- functional(ma, y~u, value=function(x) -.4*x^2)
##   d <- sim(ms,50)#,seed=1)
##   ea <- estimate(ma,d)
##   mb <- lvm()
##   mb1 <- nonlinear(mb,type="linear",y~u)
##   mb2 <- nonlinear(mb,type="quadratic",y~u)
##   mb3 <- nonlinear(mb,type="spline",knots=c(-3,-1,0,1,3),y~u)
##   mb4 <- nonlinear(mb,type="spline",knots=c(-3,-2,-1,0,1,2,3),y~u)
##   ff <- lapply(list(mb1,mb2,mb3,mb4),
##                function(m) function(data,...) twostage(ma,m,data=data,st.derr=FALSE))
##   a <- lava:::cv(ff,data=d,rep=1,K=2)
## }
