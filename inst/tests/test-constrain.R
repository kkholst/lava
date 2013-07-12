context("Constraints")

test_that("constrain (Fishers z-transform)",{
  m <- lvm(c(y1[m1:v1],y2[m2:v2])~x)
  covariance(m,y1~y2) <- "C"
  d <- sim(m,100)
  e <- estimate(m,d)
  constrain(e,rho~C+v1+v2) <-
    function(x) x[1]/(x[2]*x[3])^0.5
  cc1 <- correlation(e)
  cc2 <- constraints(e)
  expect_equivalent(cc2["rho",1],cc1["y1~y2",1])
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
  expect_equal(cc2["z",2],0.1)
  expect_equivalent(cc2["inv(z)",1],cc1["y1~y2",1])
})


test_that("Multiple group constraints", {
  data(twindata)
  twinwide <- reshape(twindata,direction="wide",
                      idvar="id",timevar="twinnum")
  l <- lvm(~bw.1+bw.2)
  covariance(l) <- bw.1 ~ bw.2
  e <- estimate(l,subset(twinwide,zyg.1=="MZ"))
  B <- cbind(1,-1); colnames(B) <- c("bw.1,bw.1","bw.2,bw.2")
  lava::compare(e,contrast=B)
  B2 <- rbind(c(1,-1,0,0),c(0,0,1,-1))
  colnames(B2) <- c("bw.1","bw.2","bw.1,bw.1","bw.2,bw.2")
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
})

## test_that("text",{
##   ##  expect_output(g,"p=12")
## })


