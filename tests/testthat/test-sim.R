context("Simulation")

test_that("Missing", {
    m <- lvm(y~1)
    m <- Missing(m,y~1,r~x)
    set.seed(1)
    d <- simulate(m,1e3,seed=1)
    expect_equal(sum(d$r),sum(!is.na(d$y0)))

    g <- glm(r~x,data=d,family=binomial)
    expect_true(all.equal(coef(g),c(0,1),tolerance=0.2,check.attributes=FALSE))
})


test_that("sim.default", {
    m <- lvm(y~x+e)
    distribution(m,~y) <- 0
    distribution(m,~x) <- uniform.lvm(a=-1.1,b=1.1)
    transform(m,e~x) <- function(x) (1*x^4)*rnorm(length(x),sd=1)

    onerun <- function(iter=NULL,...,n=2e3,b0=1,idx=2) {
        d <- sim(m,n,p=c("y~x"=b0))
        l <- lm(y~x,d)
        res <- c(coef(summary(l))[idx,1:2],
                 confint(l)[idx,],
                 estimate(l,only.coef=TRUE)[idx,2:4])
        names(res) <- c("Estimate","Model.se","Model.lo","Model.hi",
                        "Sandwich.se","Sandwich.lo","Sandwich.hi")
        res
    }

    val <- sim(onerun,R=2,b0=1,n=10,messages=0,mc.cores=1)
    expect_true(nrow(val)==2)
    val <- sim(val,R=2,b0=1,n=10,mc.cores=1,messages=0) ## append results
    expect_true(nrow(val)==4)

    s1 <- summary(val,estimate=c(1,1),confint=c(3,4,6,7),true=c(1,1),names=c("Model","Sandwich"))
    expect_true(length(grep("Coverage",rownames(s1)))>0)
    expect_equivalent(colnames(s1),c("Model","Sandwich"))

})


