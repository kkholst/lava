context("Inference")

test_that("estimate.default", {
    m <- lvm(c(y1,y2)~x+z,y1~~y2)
##    set.seed(1)
    d <- sim(m,20)

    l1 <- lm(y1~x+z,d)
    l2 <- lm(y2~x+z,d)
    ll <- merge(l1,l2)
    testthat::expect_equivalent(ll$coefmat[,1],c(coef(l1),coef(l2)))

    e1 <- estimate(l1)
    f1 <- estimate(l1,function(x) x^2, use=2)
    testthat::expect_true(coef(l1)["x"]^2==f1$coefmat[1])

    e1b <- estimate(NULL,coef=coef(l1),vcov=vcov(estimate(l1)))
    e1c <- estimate(NULL,coef=coef(l1), IC=IC(l1))
    testthat::expect_equivalent(vcov(e1b),vcov(e1c))
    testthat::expect_equivalent(var_ic(IC(e1)), vcov(e1b))

    f1b <- estimate(e1b,function(x) x^2)
    testthat::expect_equivalent(f1b$coefmat[2,,drop=FALSE],f1$coefmat)

    h1 <- estimate(l1,cbind(0,1,0))
    testthat::expect_true(h1$coefmat[,5]==e1$coefmat["x",5])

    ## GEE
    if (requireNamespace("geepack",quietly=TRUE)) {
        dd <- reshape(d, direction='long', varying=list(c('y1','y2')), v.names='y')
        dd <- dd[order(dd$id),]
        ## dd <- mets::fast.reshape(d)
        l <- lm(y~x+z,dd)
        g1 <- estimate(l,id=dd$id)
        g2 <- geepack::geeglm(y~x+z,id=dd$id,data=dd)
        testthat::expect_equivalent(g1$coefmat[,c(1,2,5)],
                          as.matrix(summary(g2)$coef[,c(1,2,4)]))
    }

    ## Several parameters
    e1d <- estimate(l1, function(x) list("X"=x[2],"Z"=x[3]))
    testthat::expect_equivalent(e1d$coefmat,e1$coefmat[-1,])

    testthat::expect_true(rownames(estimate(l1, function(x) list("X"=x[2],"Z"=x[3]),keep="X")$coefmat)=="X")
    testthat::expect_true(rownames(estimate(l1, labels=c("a"), function(x) list("X"=x[2],"Z"=x[3]),keep="X")$coefmat)=="a")


    a0 <- estimate(l1,function(p,data) p[1]+p[2]*data[,"x"], average=TRUE)
    a1 <- estimate(l1,function(p,data) p[1]+p[2]*data[,"x"]+p[3], average=TRUE)
    a <- merge(a0,a1,labels=c("a0","a1"))
    estimate(a,diff)
    testthat::expect_equivalent(estimate(a,diff)$coefmat,e1$coefmat[3,,drop=FALSE])


    stack

})
