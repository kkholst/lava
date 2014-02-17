context("Inference")

test_that("Effects",{
    m <- lvm()
    regression(m) <- c(y1,y2,y3)~u; latent(m) <- ~u
    regression(m) <- c(z1,z2,z3)~v; latent(m) <- ~v
    regression(m) <- u~v
    regression(m) <- c(u,v,z3,y1)~x
    d <- sim(m,100)
    start <- c(rep(0,6),rep(1,17))
    suppressWarnings(e <- estimate(m,d,control=list(iter.max=0,start=start)))
    f <- coef(effects(e,y1~x))
    expect_true(all(f[,2]>0)) ## Std.err
    expect_equal(f["Total",1],3) 
    expect_equal(f["Direct",1],1)
    f2 <- coef(effects(e,u~v))
    expect_equal(f2["Total",1],1)
    expect_equal(f2["Direct",1],1)
    expect_equal(f2["Indirect",1],0)
})

test_that("Profile confidence limits", {
    m <- lvm(y~b*x)
    constrain(m,b~psi) <- identity
    set.seed(1)
    d <- sim(m,100)
    e <- estimate(m, d)
    ci0 <- confint(e,3)
    ci <- confint(e,3,profile=TRUE)
    expect_true(mean((ci0-ci)^2)<0.1)
}
