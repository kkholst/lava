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
})

test_that("IV-estimator", {
    m <- lvm(c(y1,y2,y3)~u); latent(m) <- ~u    
    set.seed(1)
    d <- sim(m,100)
    e0 <- estimate(m,d)
    e <- estimate(m,d,estimator="iv") ## := MLE
    expect_true(mean((coef(e)-coef(e0))^2)<1e-9)
})

test_that("glm-estimator", {         
    m <- lvm(y~x+z)
    regression(m) <- x~z
    distribution(m,~y+z) <- binomial.lvm("logit")
    set.seed(1)
    d <- sim(m,1e3)
    head(d)
    e <- estimate(m,d,estimator="glm")
    c1 <- coef(e,2)[c("y","y~x","y~z"),1:2]
    c2 <- estimate(glm(y~x+z,d,family=binomial))$coefmat[,1:2]  
    expect_equivalent(c1,c2)
})


test_that("gaussian", {
    m <- lvm(y~x)
    d <- simulate(m,100,seed=1)
    S <- cov(d[,vars(m),drop=FALSE])
    mu <- colMeans(d[,vars(m),drop=FALSE])
    f <- function(p) lava:::gaussian_objective.lvm(p,x=m,S=S,mu=mu,n=nrow(d))
    g <- function(p) lava:::gaussian_score.lvm(p,x=m,n=nrow(d),data=d,indiv=TRUE)
    s1 <- numDeriv::grad(f,c(0,1,1))
    s2 <- g(c(0,1,1))
    expect_equal(s1,-colSums(s2),tolerance=0.1)
})


test_that("Association measures", {
    P <- matrix(c(0.25,0.25,0.25,0.25),2)
    a1 <- lava:::assoc(P)
    expect_equivalent(-log(0.25),a1$H)
    expect_true(with(a1, all(c(kappa,gamma,MI,U.sym)==0)))
    
    p <- lava:::prob.normal(sigma=diag(nrow=2),breaks=c(-Inf,0),breaks2=c(-Inf,0))[1]
    expect_equal(p[1],0.25)

    ## q <- qnorm(0.75)
    ## m <- ordinal(lvm(y~x),~y, K=3)#, breaks=c(-q,q))
    ## normal.threshold(m,p=c(0,1,2))
})


test_that("Bootstrap", {
    y <- rep(c(0,1),each=5)
    x <- 1:10
    e <- estimate(y~x)
    B1 <- bootstrap(e,R=2,silent=TRUE)
    B2 <- bootstrap(e,R=2,silent=TRUE,bollenstine=TRUE)
    expect_false(B1$bollenstine)
    expect_true(B2$bollenstine)
    expect_true(nrow(B1$coef)==2)
})
