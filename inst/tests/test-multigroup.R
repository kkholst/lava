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



test_that("Multiple group constraints", {
  data(twindata)
  twinwide <- reshape(twindata,direction="wide",
                      idvar="id",timevar="twinnum")
  l <- lvm(~bw.1+bw.2)
  covariance(l) <- bw.1 ~ bw.2
  e <- estimate(l,subset(twinwide,zyg.1=="MZ"))
  B <- cbind(1,-1); colnames(B) <- c("bw.1<->bw.1","bw.2<->bw.2")
  lava::compare(e,contrast=B)
  B2 <- rbind(c(1,-1,0,0),c(0,0,1,-1))
  colnames(B2) <- c("bw.1","bw.2","bw.1<->bw.1","bw.2<->bw.2")
  lava::compare(e,contrast=B2)

  l <- lvm(~bw.1+bw.2)
  covariance(l) <- bw.1 ~ bw.2
  intercept(l,~bw.1+bw.2) <- "m"
  covariance(l,~bw.1+bw.2) <- "s"
  covariance(l,bw.1~bw.2) <- "r1"
  l2 <- l
  covariance(l2,bw.1~bw.2) <- "r2"

  DZ <- subset(twinwide,zyg.1=="MZ")
  MZ <- subset(twinwide,zyg.1=="DZ")
  e <- estimate(l,MZ)
  e2 <- estimate(l2,DZ)

  parameter(l) <- ~r2
  parameter(l2) <- ~r1
  ee <- estimate(list(MZ=l,DZ=l2),list(MZ,DZ))
  constrain(ee,h~r1+r2) <- function(x) 2*(x[1]-x[2])
  ce <- constraints(ee)
  expect_true(length(coef(ee))==4)
  expect_true(nrow(ce)==1)
  expect_true(all(!is.na(ce)))
  expect_true(sum(score(ee)^2)<1e-4)
}



