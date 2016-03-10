context("Utility functions")

test_that("By", {
    require('data.table')
    By(datasets::CO2,~Treatment+Type,colMeans,~conc)
    t1 <- as.data.frame(data.table(datasets::CO2)[,mean(uptake),by=.(Treatment,Type,conc>500)])
    t2 <- By(datasets::CO2,~Treatment+Type+I(conc>500),colMeans,~uptake)
    expect_true(inherits(t2,"array"))
    expect_equivalent(sort(t2),sort(t1$V1))    
    
})


test_that("Expand", {
          dd <- Expand(iris, Sepal.Length=2:8, Species=c("virginica","setosa"))
          ##' summary(dd)
          ##'
          d0 <- datasets::warpbreaks[,c("wool","tension")]
          T <- table(d0)
          d1 <- Expand(T)
          expect_identical(dim(d0),dim(d1))
          expect_identical(table(d1),T)
})


test_that("trim", {
    expect_true(length(grep(" ",trim(" test ")))==0)    
    expect_true(length(grep(" ",trim(" t e s t ")))==3)
    expect_true(length(grep(" ",trim(" t e s t ",all=TRUE)))==0)
})


test_that("addattr", {
    m <- lvm(y~x)
    g <- plot(m,noplot=TRUE)
})

