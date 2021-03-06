context("Constraints")

test_that("Simple linear constraint",{
    m1 <- lvm(y[m:v] ~ f(x,beta)+f(z,beta2))
    constrain(m1,beta2~psi) <- function(x) 2*x
    testthat::expect_output(summary(m1),"Non-linear constraints:")

    lava:::matrices.lvm(m1,1:2,0,3)
    d1 <- sim(m1,100)
    e1 <- estimate(m1,d1)
    testthat::expect_output(print(e1),"y~x")
    testthat::expect_warning(e1,NA)

    testthat::expect_true((constraints(e1)[1]-coef(lm(y~x+z,d1))["z"])^2<1e-9)
    s <- summary(e1)
    testthat::expect_output(print(s),"Non-linear constraints:")

    testthat::expect_equivalent(dim(coef(s)), c(length(coef(e1))+1,4))
    testthat::expect_equivalent(dim(coef(e1,2)), c(length(coef(e1)),4))
})


test_that("constrain (Fishers z-transform)",{
    set.seed(1)
    m <- lvm(c(y1[m1:v1],y2[m2:v2])~x)
    covariance(m,y1~y2) <- "C"
    d <- sim(m,100)
    e <- estimate(m,d)
    constrain(e,rho~C+v1+v2) <-
        function(x) x[1]/(x[2]*x[3])^0.5
    cc1 <- coef(summary(correlation(e)))
    cc2 <- constraints(e)
    testthat::expect_equivalent(cc2["rho",1],cc1["y1~y2",1])
    constrain(e,z~C+v1+v2) <- function(x) {
        f <- function(p) p[1]/sqrt(p[2]*p[3])
        res <- atanh(f(x))
        df <- function(p) c(1/sqrt(p[2]*p[3]), -f(p)/(2*p[2]), -f(p)/(2*p[3]))
        datanh <- function(r) 1/(1-r^2)
        attributes(res)$grad <- function(p) datanh(f(p))*df(p)
        attributes(res)$inv <- tanh
        return(res)
    }
    cc2 <- constraints(e)
    testthat::expect_equal(cc2["z",2],0.1)
    testthat::expect_equivalent(cc2["inv(z)",1],cc1["y1~y2",1])
})


test_that("Non-linear in exogenous variables", {
    d <- data.frame(x=1:5,y=c(0.5,1.5,2,3,3.5))
    m <- lvm(y[m] ~ 1)
    addvar(m) <- ~x
    parameter(m) <- ~a+b
    constrain(m,m~a+b+x) <- function(z) z[1]+z[2]*z[3]
    e <- estimate(m,d,control=list(method="NR"))
    testthat::expect_true(mean(coef(lm(y~x,d))-coef(e)[c("a","b")])^2<1e-3)
})


if (lava:::versioncheck('mets', c(1,0)))
test_that("Probit constraints", {
    x <- transform(data.frame(lava:::rmvn0(1000,sigma=0.5*diag(2)+0.5)),
                   X1=as.numeric(cut(X1,breaks=3))-1,X2=as.numeric(cut(X2,breaks=3))-1)
    m <- covariance(lvm(),X1~X2)
    ordinal(m,K=3,constrain=list("t1","t2")) <- ~X1
    ordinal(m,K=3,constrain=list("t1","t2")) <- ~X2
    ##        e <- estimate(m,x)
    e <- estimate(list(m,m),list(x[1:500,],x[501:1000,]),estimator="normal")
    res <- estimate(e)
    testthat::expect_true(length(coef(res))==4)
})


test_that("Multiple group constraints I", {
    m1 <- lvm(y[m:v] ~ f(x,beta)+f(z,beta2))
    d1 <- sim(m1,500,seed=1); d2 <- sim(m1,500,seed=2)
    ##coef(estimate(m1,d1))
    constrain(m1,beta2~psi) <- function(x) 2*x
    m2 <- lvm(y[m:v] ~ f(x,beta2) + z)
    constrain(m2,beta2~psi) <- function(x) 2*x
    mg <- multigroup(list(m1,m2),list(d1,d2))
    ee <- estimate(mg)
    testthat::expect_true(length(coef(ee))==5)
    testthat::expect_equivalent(constraints(ee)[1],2*coef(ee)["psi@1"]) # Est
    testthat::expect_equivalent(constraints(ee)[2],2*coef(ee,2)[[1]]["psi",2]) # Std.Err
})

test_that("Multiple group constraints II", {
  data("twindata",package="lava")
  twinwide <- reshape(twindata,direction="wide",
                      idvar="id",timevar="twinnum")
  l <- lvm(~bw.1+bw.2)
  covariance(l) <- bw.1 ~ bw.2
  e <- estimate(l,subset(twinwide,zyg.1=="MZ"),control=list(method="NR"))
  B <- cbind(1,-1); colnames(B) <- c("bw.1,bw.1","bw.2,bw.2")
  colnames(B) <- gsub(",",lava.options()$symbols[2],colnames(B))
  lava::compare(e,contrast=B)
  B2 <- rbind(c(1,-1,0,0),c(0,0,1,-1))
  colnames(B2) <- c("bw.1","bw.2","bw.1,bw.1","bw.2,bw.2")
  colnames(B2) <- gsub(",",lava.options()$symbols[2],colnames(B2))

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
  ## e <- estimate(l,MZ)
  ## e2 <- estimate(l2,DZ)

  parameter(l) <- ~r2
  parameter(l2) <- ~r1
  ee <- estimate(list(MZ=l,DZ=l2),list(MZ,DZ),control=list(method="NR",tol=1e-9,constrain=FALSE))
  testthat::expect_true(mean(score(ee)^2)<1e-9)

  constrain(ee,h~r2+r1) <- function(x) 2*(x[1]-x[2])
  ce <- constraints(ee)
  testthat::expect_equivalent(constraints(ee)[1],2*diff(coef(ee)[3:4]))
  testthat::expect_true(length(coef(ee))==4)
  testthat::expect_true(nrow(ce)==1)
  testthat::expect_true(all(!is.na(ce)))
})
