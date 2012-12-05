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


test_that("text",{
  ##  expect_output(g,"p=12")
})


