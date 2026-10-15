# Optimization

test_that("Optimization", {
    m <- lvm(y~x+z)
    d <- simulate(m,20,seed=1)
    e1 <- estimate(m,d,control=list(method="nlminb0"))
    e2 <- estimate(m,d,control=list(method="NR"))
    testthat::expect_equal(round(coef(e1),3),round(coef(e2),3), ignore_attr = TRUE)

    f <- function(x) x^2*log(x) # x>0
    df <- function(x) 2*x*log(x) + x
    df2 <- function(x) 2*log(x) + 3
    op <- NR(5,f,df,df2,control=list(tol=1e-40)) ## Find root
    testthat::expect_equal(round(op$par,digits=7),.6065307, ignore_attr = TRUE)
    op2 <- lava:::estimatingfunction0(5,gradient=df)
    op3 <- lava:::estimatingfunction(5,gradient=df,hessian=df2,control=list(tol=1e-40))
    testthat::expect_equal(op$par,op2$par, ignore_attr = TRUE)
    testthat::expect_equal(op$par,op3$par, ignore_attr = TRUE)
})
