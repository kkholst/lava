context("Model specification")

test_that("Linear constraints", {
    m <- lvm(c(y[m:v]~b*x))
    constrain(m,b~a) <- base::identity
}) 
