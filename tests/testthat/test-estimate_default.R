context("estimate.default")

test_that("estimate.default misc", {
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

})

# Helper function to manually compute Wald statistic
compute_wald <- function(B, p, S, null) {
  z <- (B %*% p - null)
  V <- B %*% S %*% t(B)
  q <- t(z) %*% Inverse(V) %*% z
  return(structure(q[1], df=qr(V)$rank))
}

set.seed(1)
a1 <- estimate(coef = 1,   IC = rnorm(10), id = 1:10, labels = "a1")
a2 <- estimate(coef = 2,   IC = rnorm(10), id = 1:10, labels = "a2")
a3 <- estimate(coef = 3,   IC = rnorm(10), id = 1:10, labels = "a3")
a4 <- estimate(coef = 4,   IC = rnorm(10), id = 1:10, labels = "a4")
a  <- merge(a1, a2)           # 2-dimensional
a3d <- merge(a1, a2, a3)      # 3-dimensional
a4d <- merge(a1, a2, a3, a4)  # 4-dimensional


test_that("summary.estimate compared with estimate", {
  B <- rbind(c(1,-1, 0), c(0, 1,-1), c(1,0,-1))
  null <- c(1,2,3)
  q <- compute_wald(B, coef(a3d), vcov(a3d), null)
  df <- attr(q, "df")
  e1 <- estimate(a3d, f=B)
  e1b <- estimate(a3d, f=B, null=null)
  e1a <- estimate(e1, f=diag(3), null=null)
  expect_equal(unname(df), unname(e1a$compare$parameter))
  expect_equal(unname(df), unname(e1b$compare$parameter))
  expect_true(abs(pchisq(q, df=df, lower.tail=FALSE) - e1a$compare$p.value)<1e-16)
  expect_true(abs(pchisq(q, df=df, lower.tail=FALSE) - e1b$compare$p.value)<1e-16)
  expect_true(abs(q - e1a$compare$statistic) < 1e-9)
  expect_true(abs(q - e1b$compare$statistic) < 1e-9)
  # function spec.
  e2 <- estimate(a3d, function(p) c(p[1]-p[2], p[2]-p[3], p[1]-p[3]))
  e2b <- estimate(a3d, function(p) c(p[1]-p[2], p[2]-p[3], p[1]-p[3]), null=null)
  e2a <- estimate(e2, f=diag(3), null=null)
  expect_equal(unname(df), unname(e2a$compare$parameter))
  expect_equal(unname(df), unname(e2b$compare$parameter))
  expect_true(abs(pchisq(q, df=df, lower.tail=FALSE) - e2a$compare$p.value)<1e-16)
  expect_true(abs(pchisq(q, df=df, lower.tail=FALSE) - e2b$compare$p.value)<1e-16)
  expect_true(abs(q - e2a$compare$statistic) < 1e-9)
  expect_true(abs(q - e2b$compare$statistic) < 1e-9)
  # summary
  e3 <- summary(e1, null=null)
  expect_equal(unname(df), unname(e3$compare$parameter))
  expect_true(abs(pchisq(q, df=df, lower.tail=FALSE) - e3$compare$p.value)<1e-16)
  expect_true(abs(q - e3$compare$statistic) < 1e-9)
})

test_that("summary.estimate compared with estimate", {
  B <- rbind(c(2,-1, 0), c(0, 3,-1), c(1,0,-3), c(1,0,0))
  null <- c(1,0,1,0)
  q <- compute_wald(B, coef(a3d), vcov(a3d), null)
  df <- attr(q, "df")
  e1 <- estimate(a3d, f=B)
  e1b <- estimate(a3d, f=B, null=null)
  e1a <- estimate(e1, f=diag(4), null=null)
  expect_equal(unname(df), unname(e1a$compare$parameter))
  expect_equal(unname(df), unname(e1b$compare$parameter))
  expect_true(abs(pchisq(q, df=df, lower.tail=FALSE) - e1a$compare$p.value)<1e-16)
  expect_true(abs(pchisq(q, df=df, lower.tail=FALSE) - e1b$compare$p.value)<1e-16)
  expect_true(abs(q - e1a$compare$statistic) < 1e-9)
  expect_true(abs(q - e1b$compare$statistic) < 1e-9)
  # summary
  e3 <- summary(e1, null=null)
  expect_equal(unname(df), unname(e3$compare$parameter))
  expect_true(abs(pchisq(q, df=df, lower.tail=FALSE) - e3$compare$p.value)<1e-16)
  expect_true(abs(q - e3$compare$statistic) < 1e-9)
})

test_that("1D: Identity contrast, null = 0", {
  B    <- matrix(1, nrow = 1, ncol = 1)
  null <- 0
  e    <- estimate(a1, f = B, null = null)
  p <- coef(a1)
  S <- vcov(a1)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
})

test_that("1D: Identity contrast, null = 1", {
  B    <- matrix(1, nrow = 1, ncol = 1)
  null <- 1
  e <- estimate(a1, f = B, null = null)
  p <- coef(a1)
  S <- vcov(a1)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
})

test_that("1D: Scalar contrast (scaling), null = 0", {
  B    <- matrix(2, nrow = 1, ncol = 1)
  null <- 0
  e <- estimate(a1, f = B, null = null)
  p <- coef(a1)
  S <- vcov(a1)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
})

test_that("2D: Identity contrast, null = c(0, 0)", {
  B    <- diag(2)
  null <- c(0, 0)
  e <- estimate(a, f = B, null = null)
  p <- coef(a)
  S <- vcov(a)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
})

test_that("2D: Identity contrast, null equals true coefs", {
  B    <- diag(2)
  null <- c(1, 2)
  e <- estimate(a, f = B, null = null)
  p <- coef(a)
  S <- vcov(a)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
  # Wald statistic should be very small (close to 0)
  # since null = true coefs
  expect_lt(e$compare$statistic, 1e-10)
})

test_that("2D: Identity contrast, null = c(3, 4)", {
  B    <- diag(2)
  null <- c(3, 4)
  e <- estimate(a, f = B, null = null)
  p <- coef(a)
  S <- vcov(a)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
})

test_that("2D: Difference contrast c(1,-1), null = 0", {
  B    <- matrix(c(1, -1), nrow = 1)
  null <- 0
  e <- estimate(a, f = B, null = null)
  p <- coef(a)
  S <- vcov(a)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
})

test_that("2D: Difference contrast c(1,-1), null = -1", {
  B    <- matrix(c(1, -1), nrow = 1)
  null <- -1  # Testing H0: a1 - a2 = -1 (true, since 1 - 2 = -1)
  e <- estimate(a, f = B, null = null)
  p <- coef(a)
  S <- vcov(a)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
  # Wald stat ~ 0 since null equals true difference
  expect_lt(e$compare$statistic, 1e-10)
})

test_that("2D: Sum contrast c(1,1), null = 3", {
  B    <- matrix(c(1, 1), nrow = 1)
  null <- 3  # True sum = 1 + 2 = 3
  e <- estimate(a, f = B, null = null)
  p <- coef(a)
  S <- vcov(a)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
  # Wald stat ~ 0
  expect_lt(e$compare$statistic, 1e-10)
})

test_that("2D: Scaling contrast, null = c(2, 6)", {
  B    <- 2 * diag(2)
  null <- c(2, 6)  # True: 2*c(1,2) = c(2,4), so null != true
  e <- estimate(a, f = B, null = null)
  p <- coef(a)
  S <- vcov(a)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
})

test_that("3D: Identity contrast, null = true coefs", {
  B    <- diag(3)
  null <- c(1, 2, 3)
  e <- estimate(a3d, f = B, null = null)
  p <- coef(a3d)
  S <- vcov(a3d)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
  expect_lt(e$compare$statistic, 1e-10)
})

test_that("3D: Identity contrast, null = c(0, 0, 0)", {
  B    <- diag(3)
  null <- c(0, 0, 0)
  e <- estimate(a3d, f = B, null = null)
  p <- coef(a3d)
  S <- vcov(a3d)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
})

test_that("3D: Pairwise differences, null = c(-1, -2)", {
  # H0: a1-a2 = -1, a2-a3 = -1 (true: 1-2=-1, 2-3=-1)
  B <- matrix(c(
    1, -1,  0,
    0,  1, -1
  ), nrow = 2, byrow = TRUE)
  null <- c(-1, -1)
  e <- estimate(a3d, f = B, null = null)
  p <- coef(a3d)
  S <- vcov(a3d)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
  expect_lt(e$compare$statistic, 1e-10)
})

test_that("3D: Single-row contrast, null = 6", {
  # H0: a1 + a2 + a3 = 6 (true: 1+2+3=6)
  B    <- matrix(c(1, 1, 1), nrow = 1)
  null <- 6
  e <- estimate(a3d, f = B, null = null)
  p <- coef(a3d)
  S <- vcov(a3d)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
  expect_lt(e$compare$statistic, 1e-10)
})

test_that("4D: Identity contrast, null equals true coefs", {
  B    <- diag(4)
  null <- c(1, 2, 3, 4)
  e <- estimate(a4d, f = B, null = null)
  p <- coef(a4d)
  S <- vcov(a4d)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
  expect_lt(e$compare$statistic, 1e-10)
})

test_that("4D: Non-square contrast (2x4), null = c(0, 0)", {
  # H0: a1 - a2 = 0, a3 - a4 = 0
  B <- matrix(c(
    1, -1,  0,  0,
    0,  0,  1, -1
  ), nrow = 2, byrow = TRUE)
  null <- c(0, 0)
  e <- estimate(a4d, f = B, null = null)
  p <- coef(a4d)
  S <- vcov(a4d)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
})

test_that("4D: Non-square contrast, null equals true values", {
  # H0: a1-a2 = -1, a3-a4 = -1 (true: 1-2=-1, 3-4=-1)
  B <- matrix(c(
    1, -1,  0,  0,
    0,  0,  1, -1
  ), nrow = 2, byrow = TRUE)
  null <- c(-1, -1)
  e <- estimate(a4d, f = B, null = null)
  p <- coef(a4d)
  S <- vcov(a4d)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
  expect_lt(e$compare$statistic, 1e-10)
})

test_that("Degrees of freedom equals nrow(B)", {
  B    <- diag(2)
  null <- c(0, 0)
  e    <- estimate(a, f = B, null = null)
  expect_identical(unname(e$compare$parameter), nrow(B))
  B    <- matrix(c(1, -1, 0, 0, 1, -1), nrow = 2, byrow = TRUE)
  null <- c(0, 0)
  e    <- estimate(a3d, f = B, null = null)
  expect_equal(unname(e$compare$parameter), nrow(B))
})
