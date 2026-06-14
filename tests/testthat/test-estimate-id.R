context("estimate.default id/cluster")

test_that("index.estimate", {
    m <- lvm(c(y1,y2)~x+z, y1~~y2)
    set.seed(1)
    d <- sim(m,20)

    l1 <- lm(y1~x+z, data=d)
    e1 <- estimate(l1)
    testthat::expect_equivalent(index(e1), rownames(d))
    testthat::expect_true(inherits(index(e1), "character"))
    V <- vcov(e1)
    index(e1) <- as.numeric(index(e1))
    expect_equivalent(vcov(e1), V)
    testthat::expect_true(inherits(index(e1), "numeric"))
})
