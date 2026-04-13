context("Graphics functions")
library("vdiffr")

test_that("color, devcoords", {
  cur <- palette()
  old <- lava:::mypal()
  testthat::expect_equivalent(col2rgb(cur),col2rgb(old))
  testthat::expect_equivalent(col2rgb(palette()),col2rgb(lava:::mypal(set=FALSE)))

  testthat::expect_equivalent(Col("red",0.5),rgb(1,0,0,0.5))
  testthat::expect_equivalent(Col(c("red","blue"),0.5),rgb(c(1,0),c(0,0),c(0,1),0.5))
  testthat::expect_equivalent(Col(c("red","blue"),c(0.2,0.5)),rgb(c(1,0),c(0,0),c(0,1),c(0.2,0.5)))
  testthat::expect_equivalent(Col(rgb(1,0,0),0.5),rgb(1,0,0,0.5))

  plot(0,xlim=c(0,1),ylim=c(0,1),type="n",ann=FALSE,axes=FALSE)
  devc1 <- devcoords()
  par(mar=c(0,0,0,0))
  plot(0,xlim=c(0,1),ylim=c(0,1),type="n",ann=FALSE,axes=FALSE)
  devc2 <- devcoords()
  figx <- c("fig.x1","fig.x2","fig.y1","fig.y2")
  devx <- c("dev.x1","dev.x2","dev.y1","dev.y2")
  testthat::expect_equivalent(devc1[figx],devc2[devx])
})

test_that("plotConf", {
  skip_on_cran()

  set.seed(1)
  x <- rnorm(50)
  y <- rnorm(50,x)
  z <- rbinom(50,1,0.5)
  d <- data.frame(y,z,x)
  l <- lm(y~x*z)
  myplot <- function(...) {
    plotConf(l,var1="x",var2="z",col=c("black","blue"),alpha=0.5,legend=FALSE)
  }
  expect_doppelganger("plotConf1", myplot)
  ##par(mar=c(0,0,0,0))
  ## newd <- data.frame(x=seq(min(x),max(x),length.out=100))
  ## l0 <- lm(y~x,subset(d,z==0))
  ## ci0 <- predict(l0,newdata=newd,interval="confidence")
  ## l1 <- lm(y~x,subset(d,z==1))
  ## ci1 <- predict(l1,newdata=newd,interval="confidence")
  ## plot(y~x,col=c("black","blue")[z+1],pch=16,ylim=c(min(ci0,ci1,y),max(ci0,ci1,y)))
  ## lines(newd$x,ci0[,1],col="black",lwd=2)
  ## lines(newd$x,ci1[,1],col="blue",lwd=2)
  ## confband(newd$x,lower=ci0[,2],upper=ci0[,3],polygon=TRUE,col=Col("black",0.5),border=FALSE)
  ## confband(newd$x,lower=ci1[,2],upper=ci1[,3],polygon=TRUE,col=Col("blue",0.5),border=FALSE)
  ## points(y~x,col=c("black","blue")[z+1],pch=16)
  #testthat::expect_true(grcompare(d1,d2,threshold=5))
  l <- lm(y~z)
  myplot <- function(...) {
    plotConf(l,var2="z",var1=NULL,jitter=0,col="black",alpha=0.5,xlim=c(.5,2.5),ylim=range(y))
  }
  expect_doppelganger("plotConf2", myplot)
  ## par(mar=c(0,0,0,0))
  ## plot(y~I(z+1),ylim=range(y),xlim=c(0.5,2.5),pch=16,col=Col("black",0.5))
  ## l0 <- lm(y~-1+factor(z))
  ## confband(1:2,lower=confint(l0)[,1],upper=confint(l0)[,2],lwd=3,
  ##          center=coef(l0))
  })

test_that("forestplot", {
  skip_on_cran()

  set.seed(1)
  K <- 20
  est <- rnorm(K); est[c(3:4,10:12)] <- NA
  se <- runif(K,0.2,0.4)
  x <- cbind(est,est-2*se,est+2*se,runif(K,0.5,2))
  rownames(x) <- unlist(lapply(letters[seq(K)],function(x) paste(rep(x,4),collapse="")))
  rownames(x)[which(is.na(est))] <- ""
  signif <- sign(x[,2])==sign(x[,3])
  myplot <- function(...) forestplot(x)
  expect_doppelganger("forestplot1", myplot)
})

test_that("plot.sim", {
  skip_on_cran()
  skip_on_ci() # Skips the test on GitHub Actions

  set.seed(1)
  val <- as.sim(data.frame(a=rnorm(10), a.se=1, b=rnorm(10), b.se=1))
  myplot <- function(...) {
    function() plot(val, ...)
  }
  expect_doppelganger(
    "plot.sim-1",
    myplot(estimate=c(1,3), se=c(2,4), alpha=0.2,
           legend=NULL, angle=c(45,-45))
  )
  ## expect_doppelganger(
  ##   "plot.sim-2",
  ##   myplot(estimate=1, se=2, true=0, legend=NULL,
  ##          col=c("darkred", "gray"), alpha=1,
  ##          density=c(2,2), angle=c(45,-45))
  ## )
  ## expect_doppelganger(
  ##   "plot.sim-3",
  ##   myplot(estimate=c(1,3), se=c(2,4), true=0,
  ##          plot.type = "single", alpha=1, legend=NULL,
  ##          col=1:2,
  ##          density=c(2,2), angle=c(45,-45))
  ## )
  mydensity <- function(...) {
    plot(0, xlim=c(-3,3), ylim=c(0,1), type="n")
    density(val, estimate=1, legend=NULL, add=TRUE)
  }
  expect_doppelganger("density.sim-1", mydensity)
})

test_that("plot.estimate", {
  skip_on_cran()

  set.seed(1)
  e1 <- estimate(coef=1, IC=1:10, id=1:10)
  e2 <- estimate(coef=1.1, IC=0:9, id=1:10)
  myplot <- function(...) plot(c(e1, e2), null=1)
  val <- expect_doppelganger("plot.estimate-1", myplot)
})

test_that("spaghetti", {
  skip_on_cran()

  K <- 5
  y <- "y"%++%seq(K)
  m <- lvm()
  regression(m,y=y,x=~u) <- 1
  regression(m,y=y,x=~s) <- seq(K)-1
  regression(m,y=y,x=~x) <- "b"
  d <- sim(m,5,seed=1)
  dd <- mets::fast.reshape(d);
  dd$num <- dd$num+rnorm(nrow(dd),sd=0.5) ## Unbalanced
  myplot <- function(...) {
    spaghetti(y~num,dd,id="id",lty=1,col=Col(1,.4),trend=TRUE,trend.col="darkblue")
  }
  expect_doppelganger("spaghetti-1", myplot)
})

test_that("plot.lvm", {
  skip_on_cran()

  ## TODO
  m <- lvm(y~1*u[0:1],u~1*x)
  latent(m) <- ~u
  ## expect_doppelganger("plot.lvm-1", {
  ##   plot(m)
  ## })
  ## d <- sim(m,20,seed=1)
  ## e <- estimate(m,d)
  ## expect_doppelganger("plot.lvmfit-1", {
  ## plot(e)
  ## })
  ## expect_doppelganger("plot.lvmfit-2", {
  ##   plot(lava:::beautify(m))
  ## })
  g <- plot(lava:::beautify(m), noplot=TRUE)
  testthat::expect_true(inherits(g, "graphNEL"))
  g <- igraph.lvm(m)
  testthat::expect_true(inherits(g,"igraph"))
})

test_that("curly", {
  skip_on_cran()

  myplot <- function() {
  plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
  curly(x=c(1,0),y=c(0,1),lwd=2,text="a")
  curly(x=c(1,0),y=c(0,1),lwd=2,text="b",theta=pi)
  curly(x=-0.5,y=0,shape=1,theta=pi,text="c")
  curly(x=0,y=0,shape=1,theta=0,text="d")
  curly(x=0.5,y=0,len=0.2,theta=pi/2,col="blue",lty=2)
  curly(x=0.5,y=-0.5,len=0.2,theta=-pi/2,col="red",shape=1e3,text="e")
  }
  expect_doppelganger("curly-1", myplot)
})


test_that("ksmooth", {
  skip_on_cran()

  ## TODO
  expect_true(TRUE)
})

test_that("images", {
  skip_on_cran()

  ## TODO
  expect_true(TRUE)
})

test_that("labels,edgelabels", {
  skip_on_cran()

  ## TODO
  expect_true(TRUE)
})

test_that("colorbar", {
  skip_on_cran()

  ## TODO
  expect_true(TRUE)
})

test_that("fplot", {
  skip_on_cran()

  ## TODO
  expect_true(TRUE)
})

test_that("interactive", {
  skip_on_cran()

  ## TODO
  expect_true(TRUE)
})

test_that("pdfconvert", {
  skip_on_cran()

  ## TODO
  expect_true(TRUE)
})

test_that("logo", {
  skip_on_cran()

  myplot <- function(...) lava:::lava(w=10, seed=42)
  val <- expect_doppelganger("logo", myplot)
})
