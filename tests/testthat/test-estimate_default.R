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
# Generate mean-zero random IC matrix (for internal testing)
center_ic <- function(n, p = 1) {
  x <- matrix(rnorm(n * p), n, p)
  scale(x, center = TRUE, scale = FALSE)
}
a1 <- estimate(coef = 1,   IC = center_ic(10), id = 1:10, labels = "a1")
a2 <- estimate(coef = 2,   IC = center_ic(10), id = 1:10, labels = "a2")
a3 <- estimate(coef = 3,   IC = center_ic(10), id = 1:10, labels = "a3")
a4 <- estimate(coef = 4,   IC = center_ic(10), id = 1:10, labels = "a4")
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

test_that("coef.estimate warns when back.transform is set and returns untransformed coefs", {
  e_trans <- estimate(a1, back.transform = exp)
  expect_warning(coef(e_trans), "back.transform")
  suppressWarnings(expect_equal(coef(e_trans), coef(a1)))

  # verify that estimate.summary casts only a single argument when used
  # for a back.transformed object
  warns <- character(0)
  res <- withCallingHandlers(
    summary(e_trans),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(sum(grepl("back.transform", warns)), 1L)
  expect_equal(res$coefmat, summary(a1)$coefmat)
})

test_that("estimate warns when user-supplied IC has non-zero mean", {
  set.seed(42)
  ic_bad <- rnorm(50, mean = 10)
  expect_warning(
    estimate(coef = c(a = 1), IC = ic_bad, id = 1:50),
    "mean zero"
  )
  ## IC with empirical mean zero — no warning
  ic_good <- center_ic(50)
  expect_no_warning(
    estimate(coef = c(a = 1), IC = ic_good, id = 1:50)
  )
})

test_that("merge.estimate warns when input IC has non-zero mean", {
  set.seed(43)
  ic_bad <- rnorm(50, mean = 5)
  a <- suppressWarnings(estimate(coef = c(a = 1), IC = ic_bad, id = 1:50))
  b <- estimate(coef = c(b = 2), IC = center_ic(50), id = 1:50)
  expect_warning(merge(a, b), "mean zero")
})

test_that("IC mean-zero warning can be suppressed via lava.options", {
  set.seed(44)
  ic_bad <- rnorm(50, mean = 10)
  old <- lava.options(check.ic = FALSE)
  on.exit(lava.options(old))
  expect_no_warning(
    estimate(coef = c(a = 1), IC = ic_bad, id = 1:50)
  )
})

test_that("estimate.default parsedesign dispatch", {
  # Single character name selects that coefficient.
  e <- estimate(a3d, "a2")
  expect_equal(e$coefmat[, 1], 2)
  ef <- estimate(a3d, \(x) x["a2"])
  expect_equal(e$coefmat, ef$coefmat)

  # Multiple symbolic ... args: arithmetic on quoted names captured unevaluated.
  e <- estimate(a3d, "a1", "a2" - "a1", 2 * "a3" - 3 * "a1")
  expect_equal(e$coefmat[1, 1], 1)
  expect_equal(e$coefmat[2, 1], 2 - 1)
  expect_equal(e$coefmat[3, 1], 2 * 3 - 3 * 1)

  # Numeric (non-matrix) f treated as parameter index.
  e1 <- estimate(a3d, 2, 3)
  e2 <- estimate(a3d, "a2", "a3")
  expect_equal(e1$coefmat[, 1], e2$coefmat[, 1])

  # Multi-match wildcard "a*" produces pairwise contrasts (first vs rest).
  e_multi <- estimate(a3d, "a*")
  expect_equal(unname(e_multi$coefmat[, 1]), c(1 - 2, 1 - 3))

  # regex argument: pattern ".*2" diverges between glob and regex semantics.
  #   regex=FALSE: literal ".*2" matches nothing -> falls back to all coefs.
  #   regex=TRUE:  regex ".*2" matches "a2" only.
  e_glob  <- suppressWarnings(estimate(a3d, ".*2", regex = FALSE))
  e_regex <- suppressWarnings(estimate(a3d, ".*2", regex = TRUE))
  expect_equal(rownames(e_glob$coefmat), c("a1", "a2", "a3"))
  expect_equal(rownames(e_regex$coefmat), "a2")

  # Character f with user-supplied coef vector (no model object).
  pp <- c("(Intercept)" = 1.0, x = 2.0, z = -0.5)
  V  <- diag(3)
  dimnames(V) <- list(names(pp), names(pp))
  e <- estimate(f = "x", coef = pp, vcov = V)
  expect_equal(e$coefmat[, 1], 2.0)

  # Character contrasts combine with null hypothesis vector.
  e <- estimate(a3d, "a1", "a1", null = c(0, 1))
  expect_true(e$coefmat[1, "P-value"] != e$coefmat[2, "P-value"])
})

test_that("estimate.default keep with regex=TRUE", {
  e0 <- estimate(a3d, keep = ".*2", regex = TRUE)
  expect_equal(rownames(e0$coefmat), "a2")
  # the regex behavior differs from the above tests when supplying strings
  # to obtain contrasts (no matches return all coefficients)
  e1 <- estimate(a3d, keep = ".*2") # no literal matches return object with NAs
  expect_true(all(is.na(e1$coefmat)))
  expect_true(nrow(e1$coefmat) == 1)
})

test_that("robust argument backwards compatibility", {
  d <- data.frame(y = rnorm(50), x = rnorm(50))
  g <- lm(y ~ x, data = d)

  # Both robust=TRUE and robust=FALSE emit a deprecation warning
  e0 <- expect_warning(estimate(g, robust = FALSE), "deprecated and ignored")
  e1 <- expect_warning(estimate(g, robust = TRUE), "deprecated and ignored")
  expect_equal(e0$coefmat, e1$coefmat)

  e <- estimate(g)
  # The robust argument is ignored: results are identical to the default
  # (sandwich SEs)
  expect_equal(e$coefmat, e0$coefmat)

  expect_equal(e$coefmat, e0$coefmat)

  # model-based SE can be obtained either via logical variable or supplying
  # covariance matrix
  e0m <- estimate(g, vcov = TRUE)
  expect_false(all(e$coefmat == e0m$coefmat))
  e1m <- estimate(g, vcov = vcov(g))
  expect_equal(e0m$coefmat, e1m$coefmat)
})

test_that("c.summary.estimate concatenates coefficient matrices", {
  s1 <- summary(a1)
  s2 <- summary(a2)
  s3 <- summary(a3)
  cc <- c(s1, s2, s3)

  # Result is a summary.estimate object
  expect_s3_class(cc, "summary.estimate")

  # coefmat is the row-bind of the inputs (preserving rows and columns)
  expect_equal(rownames(cc$coefmat), c("a1", "a2", "a3"))
  expect_equal(colnames(cc$coefmat), colnames(s1$coefmat))
  expect_equivalent(cc$coefmat, rbind(s1$coefmat, s2$coefmat, s3$coefmat))

  # coef.summary.estimate returns the combined coefmat
  expect_equal(coef(cc), cc$coefmat)

  # original summaries are retained in the objects element
  expect_length(cc$objects, 3L)
  expect_equal(cc$objects[[2]]$coefmat, s2$coefmat)
})

test_that("c.summary.estimate with single argument returns input unchanged", {
  s1 <- summary(a1)
  expect_identical(c(s1), s1)
})

test_that("c.summary.estimate errors on non-summary.estimate input", {
  s1 <- summary(a1)
  expect_error(c(s1, a1), "only summary.estimate objects")
  expect_error(c(s1, 1), "only summary.estimate objects")
})

test_that("c.summary.estimate carries a custom print method", {
  cc <- c(summary(a1), summary(a2))
  expect_type(cc$print, "closure")
  out <- capture.output(print(cc))
  expect_true(any(grepl("Concatenated summary.estimate objects", out)))
  expect_true(any(grepl("a1", out)))
  expect_true(any(grepl("a2", out)))
})

test_that("only.coef argument is deprecated", {
  expect_warning(
    result <- estimate(a3d, only.coef = TRUE),
    regexp = "only.coef.*deprecated"
  )

  # Verify it still returns the coefficient matrix
  expect_equal(result, coef(estimate(a3d), mat = TRUE))
})
