context("Simulation")

test_that("Missing", {
    ## m <- lvm(y0~x01+x02+x03)
    ## m <- Missing(m,formula=x1~x01,Rformula=R1~0.3*x02+-0.7*x01,p=0.4)
    ## sim(m,10)
    m <- lvm(y~1)
    m <- Missing(m,y~1,r~x)
    set.seed(1)
    d <- simulate(m,1e3,seed=1)
    expect_equal(sum(d$r),sum(!is.na(d$y0)))

    g <- glm(r~x,data=d,family=binomial)
    expect_true(all.equal(coef(g),c(0,1),tolerance=0.2,check.attributes=FALSE))
})
