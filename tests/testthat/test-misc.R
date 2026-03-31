context("Utility functions")

test_that("By", {
    b1 <- By(datasets::CO2,~Treatment+Type,colMeans,~conc)
    b2 <- By(datasets::CO2,c('Treatment','Type'),colMeans,'conc')
    testthat::expect_equivalent(b1,b2)
    ## require('data.table')
    ## t1 <- as.data.frame(data.table(datasets::CO2)[,mean(uptake),by=.(Treatment,Type,conc>500)])
    d0 <- transform(datasets::CO2,conc500=conc>500)
    t1 <- by(d0[,"uptake"],d0[,c("Treatment","Type","conc500")],mean)
    t2 <- By(datasets::CO2,~Treatment+Type+I(conc>500),colMeans,~uptake)
    testthat::expect_true(inherits(t2,"array"))
    testthat::expect_equivalent(sort(t2),sort(t1))
})

test_that("Expand", {
    dd <- Expand(iris, Sepal.Length=2:8, Species=c("virginica","setosa"))
    testthat::expect_identical(levels(iris$Species),levels(dd$Species))
    testthat::expect_true(nrow(dd)==14)
    
    d0 <- datasets::warpbreaks[,c("wool","tension")]
    T <- table(d0)
    d1 <- Expand(T)
    testthat::expect_identical(dim(d0),dim(d1))
    testthat::expect_identical(table(d1),T)

    testthat::expect_identical(expand.grid(1:2,1:2),Expand(1:2,1:2))
    testthat::expect_identical(expand.grid(a=1:2,b=1:2),Expand(a=1:2,b=1:2))
})

test_that("%in%, %nin%, %in.open%, %in.closed%", {
  expect_equivalent(1:10 %ni% c(1,5,10),
                    !c(TRUE,FALSE,FALSE,FALSE,TRUE,
                      FALSE,FALSE,FALSE,FALSE,TRUE))

  expect_false(1 %in.open% c(1,4))
  expect_true(1 %in.closed% c(1,4))
  expect_error(1 %in.open% c(1,4,5))
  expect_error(1 %in.closed% c(1,4,5))


})

test_that("formulas", {
    f <- toformula(c('y1','y2'),'x'%++%1:5)
    ff <- getoutcome(f)
    testthat::expect_equivalent(trim(ff,all=TRUE),"c(y1,y2)")
    testthat::expect_true(length(attr(ff,'x'))==5)
})

test_that("trim", {
    testthat::expect_true(length(grep(" ",trim(" test ")))==0)    
    testthat::expect_true(length(gregexpr(" ",trim(" t e s t "))[[1]])==3)
    testthat::expect_true(length(grep(" ",trim(" t e s t ",all=TRUE)))==0)
})

test_that("Matrix operations:", {
    ## vec operator
    testthat::expect_equivalent(vec(diag(3)),c(1,0,0,0,1,0,0,0,1))
    testthat::expect_true(nrow(vec(diag(3),matrix=TRUE))==9)

    ## commutaion matrix
    A <- matrix(1:16 ,ncol=4)
    K <- commutation(A)
    testthat::expect_equivalent(K%*%as.vector(A),vec(t(A),matrix=TRUE))
})

test_that("blockdiag:", {
    ## Block diagonal
    A <- diag(3)+1
    B <- blockdiag(A,A,A,pad=NA)
    testthat::expect_equivalent(dim(B),c(9,9))
    testthat::expect_true(sum(is.na(B))==81-27)

    B <- matrix(1, ncol=2, nrow=2)
    res <- blockdiag(A, B, B)
    expect_equivalent(res[1:3, 1:3], A)
    expect_equivalent(res[4:5, 4:5], B)
    expect_equivalent(res[6:7, 6:7], B)
    expect_true(sum(res==0) == (7*7 - 3*3 - 2*2 - 2*2))
})

test_that("wrapvev", {
    testthat::expect_equivalent(wrapvec(5,2),c(3,4,5,1,2))
    testthat::expect_equivalent(wrapvec(seq(1:5),-1),c(5,1,2,3,4))
})

test_that("revdiag", {
    A <- revdiag(1:3)
    testthat::expect_equivalent(A,matrix(c(0,0,1,0,2,0,3,0,0),3))
    testthat::expect_equivalent(1:3,revdiag(A))
    revdiag(A) <- 4
    testthat::expect_equivalent(rep(4,3),revdiag(A))
    diag(A) <- 0
    offdiag(A) <- 5
    testthat::expect_true(sum(offdiag(A))==6*5)
})

test_that("Inv, matrix inverse", {
    A <- matrix(0,3,3)
    offdiag(A,type=3) <- 1:6
    B <- crossprod(A)
    
    testthat::expect_equivalent(solve(A),Inverse(A))
    testthat::expect_equivalent(det(B),attr(Inverse(B,chol=TRUE),"det"))
})

test_that("Grep", {
  d <- Grep(iris, "Sepal")
  expect_true(nrow(d) == nrow(iris))
  expect_equal(colnames(d), c("Sepal.Length", "Sepal.Width"))
  m <- as.matrix(iris[,1:4])
  d <- Grep(m, "Sepal")
  expect_true(nrow(d) == nrow(iris))
  expect_equal(colnames(d), c("Sepal.Length", "Sepal.Width"))

  d <- Grep(iris, "Sepal", subset=FALSE)
  expect_equivalent(colnames(d), c("index", "name"))
  expect_true(nrow(d)==2L)
})

test_that("All the rest", {
    testthat::expect_false(lava:::versioncheck(NULL))
    testthat::expect_true(lava:::versioncheck("lava",c(1,4,1)))

    op <- lava.options(debug=TRUE)
    testthat::expect_true(lava.options()$debug)
    ## lava.options(op)
     lava.options(debug=FALSE)

    A <- diag(2); colnames(A) <- c("a","b")    
    testthat::expect_output(lava:::printmany(A,A,2,rownames=c("A","B"),bothrows=FALSE),"a b")
    testthat::expect_output(lava:::printmany(A,A[1,,drop=FALSE],2,rownames=c("A","B"),bothrows=FALSE),"a b")
    testthat::expect_output(lava:::printmany(A,A,2,rownames=c("A","B"),name1="no.1",name2="no.2",
                            bothrows=TRUE),"no.1")
})

test_that("napass.0", {
  n <- 10
  d <- data.frame(y=rnorm(n), a=rbinom(n,1,0.5)*2-1)
  idx1 <- c(1,3,5)
  idx2 <- c(1,3,8,9)
  d$y[idx1] <- NA
  d$a[idx2] <- NA
  d0 <- na.pass0(d)
  expect_true(nrow(d0)==n)
  expect_true(all(d0$a[idx2]==0))
  expect_true(all(d0$a[idx1]==0))
})
