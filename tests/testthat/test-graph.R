context("Inference")

test_that("d-separation",{
    m <- lvm(x5 ~ x4+x3, x4~x3+x1, x3~x2, x2~x1)
    testthat::expect_true(dsep(m,x5~x1|x3+x4))
    testthat::expect_false(dsep(m,x5~x1|x2+x4))
    testthat::expect_true(dsep(m,x5~x1|x2+x3+x4))
    testthat::expect_false(dsep(m,~x1+x2+x3|x4))

    testthat::expect_true(setequal(ancestors(m,~x5),setdiff(vars(m),"x5")))    
    testthat::expect_true(setequal(ancestors(m,~x1),NULL))
    testthat::expect_true(setequal(descendants(m,~x5),NULL))
    testthat::expect_true(setequal(descendants(m,~x1),setdiff(vars(m),"x1")))
})

