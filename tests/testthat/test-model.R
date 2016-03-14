context("Model specification")

test_that("Linear constraints", {
    m <- lvm(c(y[m:v]~b*x))
    constrain(m,b~a) <- base::identity
}) 


test_that("Graph attributes", {
    require("graph")
    m <- lvm(y~x)
    g1 <- graph::updateGraph(plot(m,noplot=TRUE))
    col <- "blue"; v <- "y"
    g1 <- lava::addattr(g1,"fill",v,col)
    expect_match(col,graph::nodeRenderInfo(g1)$fill[v])
    nodecolor(m,v) <- "blue"
    g2 <- plot(m,noplot=TRUE)
    expect_match(col,graph::nodeRenderInfo(g2)$fill[v])
    expect_match(addattr(g2,"fill")["y"],"blue")
})


test_that("Basic model building blocks", {
    m <- lvm(y[m]~x)
    covariance(m) <- y~z
    expect_true(covariance(m)$rel["z","y"]==1)
    expect_true(regression(m)$rel["x","y"]==1)

    ## Children parent,nodes
    expect_match(children(m,~x),"y")
    expect_match(parents(m,~y),"x")    
    expect_equivalent(parents(m),vars(m))
    expect_equivalent(children(m),vars(m))

    ## Remove association
    cancel(m) <- y~z+x
    expect_true(covariance(m)$rel["z","y"]==0)
    expect_true(regression(m)$rel["x","y"]==0)

    ## Remove variable
    kill(m) <- ~x
    expect_equivalent(vars(m),c("y","z"))
    expect_true(intercept(m)["y"]=="m")  

    m <- lvm(c(y1,y2,y3)~x)
    d <- sim(m,50)
    e <- estimate(m,d)
    ## Equivalence
    ##equivalence(e,silent=TRUE)
    
}) 



test_that("Categorical variables", {
    m <- lvm()
    categorical(m,K=3,p=c(0.1,0.5)) <- ~x
    d1 <- simulate(m,10,seed=1)
    categorical(m,K=3) <- ~x
    d2 <- simulate(m,10,seed=1)
    expect_false(identical(d1,d2))
    
    regression(m,additive=FALSE,y~x) <- c(0,-5,5)
    d <- simulate(m,100,seed=1)
    l <- lm(y~factor(x),d)
    expect_true(sign(coef(l))[2]==-sign(coef(l))[3])
     
})
