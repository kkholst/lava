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

test_that("estimate standardization (average=TRUE)", {
  # check that g-computation works as expected for logistic regression
  sim1 <- function(n = 5000, seed = 1) {
    set.seed(seed)
    w1 <- rnorm(n)
    w2 <- rnorm(n)
    a  <- rbinom(n, 1, 0.5) # randomized trial
    lp <- 1 + a + w1 + 0.5 * w2
    y <- rbinom(n, 1, plogis(lp))
    data.frame(y = y, a = a, w1 = w1, w2 = w2)
  }
  df <- sim1()

  g <- glm(y ~ a * (w1 + w2), data=df, family=binomial)
  est <- lava::estimate(g, average = TRUE, function(p,data) {
    X1 <- model.matrix(g, data=transform(data, a=1))
    X0 <- model.matrix(g, data=transform(data, a=0))
    cbind(plogis(X1%*%p), plogis(X0%*%p))
  }) |> labels(c("y1", "y0"))
  est
  ## transform(est, cbind(1,-1), labels="ate")

  q1 <- predict(g, newdata=transform(df, a=1), type="response")
  q0 <- predict(g, newdata=transform(df, a=0), type="response")
  X1 <- model.matrix(g, data=transform(df, a=1))
  X0 <- model.matrix(g, data=transform(df, a=0))

  D1 <- numDeriv::grad(\(x) mean(plogis(X1 %*% x)), coef(g))
  D0 <- numDeriv::grad(\(x) mean(plogis(X0 %*% x)), coef(g))
  D1a <- apply(X1, 2, \(x) mean(x*q1*(1-q1)))
  D0a <- apply(X0, 2, \(x) mean(x*q0*(1-q0)))
  testthat::expect_equivalent(D1, D1a)
  testthat::expect_equivalent(D0, D0a)
  ic1 <- q1 - mean(q1) + apply(lava::IC(g), 1, \(x) sum(x*D1a))
  ic0 <- q0 - mean(q0) + apply(lava::IC(g), 1, \(x) sum(x*D0a))
  est2 <- c(y1=lava::estimate(coef=mean(q1), IC=ic1),
            y0=lava::estimate(coef=mean(q0), IC=ic0))
  testthat::expect_equivalent(coef(est), coef(est2))
  testthat::expect_equivalent(vcov(est), vcov(est2))

})

test_that("standardization with model estimated on a subset (id alignment)", {
  set.seed(1)
  n <- 300
  dat <- data.frame(w1 = rnorm(n),
                    a = rbinom(n, 1, 0.5),
                    z = rbinom(n, 1, 0.5))
  dat$y <- rbinom(n, 1, plogis(-0.5 + dat$w1 + dat$a))
  dat$id <- paste0("a", seq_len(n)) # lexicographic order != row order
  dat1 <- subset(dat, z == 1)
  g <- glm(y ~ w1 + a, data = dat1, family = binomial)
  f <- function(p, data) plogis(p[1] + p["w1"] * data[, "w1"] + p["a"])

  # target: E[U(W_1, A = 1)]
  # U fitted on {Z == 1}, U ~ E[Y|W, A, Z = 1] (logistic model)
  # manual influence function
  X <- cbind(1, dat$w1, 1) # intercept, w1, a
  q <- as.vector(plogis(X %*% coef(g)))
  D <- colMeans(X * q * (1 - q))
  ic2 <- rep(0, n)
  ic2[dat$z == 1] <- IC(g) %*% D * n / nrow(dat1)
  ic <- q - mean(q) + ic2 # following the order of dat

  e <- estimate(g, id = dat1$id)
  a <- estimate(e, f, data = dat, id = dat$id, average = TRUE)
  expect_equivalent(coef(a), mean(q))
  expect_equivalent(vcov(a), sum(ic^2) / n^2)
  expect_equivalent(IC(a), ic)
  expect_identical(index(a), dat$id)

  # plain glm: model identified by the rownames of its model frame, which do
  # not overlap with dat$id (independence)
  expect_message(
    a2 <- estimate(g, f, data = dat, id = dat$id, average = TRUE),
    "independence"
  )
  expect_equivalent(coef(a2), coef(a))
  expect_equivalent(
    vcov(a2),
    var_ic(q - mean(q)) + var_ic(IC(g) %*% D) # independence
  )
  expect_identical(index(a2), c(dat$id, rownames(dat1)))
  a3 <- estimate(g, f, data = dat, average = TRUE) # id defaults to rownames
  expect_equivalent(vcov(a3), vcov(a))

  # row order of 'data' does not matter
  ord <- sample(n)
  a4 <- estimate(e, f, data = dat[ord, ], id = dat$id[ord], average = TRUE)
  expect_equivalent(vcov(a4), vcov(a))
  expect_equivalent(IC(a4)[dat$id, ], IC(a)[dat$id, ])

  # conditional average (subset)
  # target: E[U(W_1, A = 1) | S = 1], for subset indicator S
  s <- dat$w1 > 0
  phat <- mean(s)
  m <- mean(q * s)
  D_s <- colMeans(X * q * (1 - q) * s)
  ic2_s <- rep(0, n)
  ic2_s[dat$z == 1] <- IC(g) %*% D_s * n / nrow(dat1)
  icc <- (q * s - m + ic2_s) / phat - m / phat^2 * (s - phat)
  ac <- estimate(e, f, data = dat, id = dat$id, subset = s, average = TRUE)
  expect_equivalent(coef(ac), m / phat)
  expect_equivalent(vcov(ac), sum(icc^2) / n^2)
  expect_equivalent(IC(ac)[,1], icc)
  expect_equivalent(ac$id, dat$id)

  # multiple parameters
  # target: E[U(W_1, A = 1)], different U
  f2 <- function(p, data) {
    list(p0 = plogis(p[1] + p["w1"] * data[, "w1"]),
         p1 = plogis(p[1] + p["w1"] * data[, "w1"] + p["a"]))
  }
  a5 <- estimate(e, f2, data = dat, id = dat$id, average = TRUE)
  expect_equivalent(coef(a5)[2], coef(a))
  expect_equivalent(vcov(a5)[2, 2], vcov(a))
  expect_equivalent(IC(a5)[,2], IC(a)[,1])

  # disjoint ids: independence between model and new data
  e_b <- estimate(g, id = paste0("b", seq_len(nrow(dat1))))
  expect_message(
    a6 <- estimate(e_b, f, data = dat, id = dat$id, average = TRUE),
    "independence"
  )
  expect_equivalent(
    vcov(a6),
    sum((q - mean(q))^2) / n^2 + var_ic(IC(g) %*% D)
  )
  expect_equal(nrow(IC(a6)), n + nrow(dat1))
  expect_identical(index(a6), c(dat$id, index(e_b)))

  # ids of estimate objects are used as is (not mapped via rownames of data):
  # default (rowname) ids of estimate(g) do not overlap with dat$id
  v_indep <- sum((q - mean(q))^2) / n^2 + var_ic(IC(g) %*% D)
  expect_message(
    a7 <- estimate(estimate(g), f, data = dat, id = dat$id, average = TRUE),
    "independence"
  )
  expect_equal(nrow(IC(a7)), n + nrow(dat1))
  expect_equivalent(vcov(a7), v_indep)

  # estimate object without index: rownames of the IF are the ids
  ic_g <- IC(g)
  e0 <- estimate(coef = coef(g), IC = ic_g)
  expect_identical(index(e0), rownames(ic_g))
  expect_message(
    a8 <- estimate(e0, f, data = dat, id = dat$id, average = TRUE),
    "independence"
  )
  expect_equal(nrow(IC(a8)), n + nrow(dat1))
  expect_equivalent(vcov(a8), v_indep)

  # estimate object without any ids: no link to data of different size
  # (no rownames in IC either)
  expect_error(
    estimate(estimate(g, id = NULL), f, data = dat, id = dat$id,
             average = TRUE),
    "Unable to link"
  )
})

test_that("standardization with partly overlapping ids", {
  set.seed(3)
  N <- 1000
  pop <- data.frame(w = rnorm(N), a = rbinom(N, 1, 0.5))
  pop$y <- rbinom(N, 1, plogis(-0.5 + pop$w + pop$a))
  pop$id <- paste0("u", seq_len(N))
  dm <- pop[401:1000, ] # model data
  dd <- pop[1:600, ]    # data for the standardization
  g <- glm(y ~ w + a, data = dm, family = binomial)
  f <- function(p, data) plogis(p[1] + p["w"] * data[, "w"] + p["a"])
  a <- estimate(estimate(g, id = dm$id), f, data = dd, id = dd$id,
                average = TRUE)
  expect_equal(nrow(IC(a)), N)
  expect_identical(index(a), pop$id)

  # manual: a_i = 1(i in data)(f_i - psi)/n_d + 1(i in model) D phi_i / n_m
  X <- cbind(1, dd$w, 1)
  q <- as.vector(plogis(X %*% coef(g)))
  D <- colMeans(X * q * (1 - q))
  ai <- rep(0, N)
  ai[1:600] <- (q - mean(q)) / 600
  ai[401:1000] <- ai[401:1000] + IC(g) %*% D / 600
  expect_equivalent(coef(a), mean(q))
  expect_equivalent(vcov(a), sum(ai^2))
  expect_equivalent(IC(a), N * ai)
})

test_that("standardization with unsorted and clustered ids (same data)", {
  set.seed(2)
  n <- 200
  d <- data.frame(x = rnorm(n))
  d$y <- rbinom(n, 1, plogis(d$x))
  d$id <- paste0("a", seq_len(n))
  d$cl <- rep(sample(paste0("c", 1:50)), each = 4)
  g <- glm(y ~ x, data = d, family = binomial)
  f <- function(p, data) plogis(p[1] + p[2] * data[, "x"])
  a0 <- estimate(g, f, average = TRUE)
  a1 <- estimate(g, f, id = d$id, average = TRUE)
  expect_equivalent(vcov(a1), vcov(a0))
  expect_equivalent(IC(a1), IC(a0))

  X <- model.matrix(g)
  q <- as.vector(plogis(X %*% coef(g)))
  ic <- q - mean(q) + IC(g) %*% colMeans(X * q * (1 - q))
  icc <- rowsum(ic, d$cl, reorder = FALSE) * 50 / n
  a2 <- estimate(g, f, id = d$cl, average = TRUE)
  expect_equivalent(vcov(a2), var_ic(icc))
  expect_identical(index(a2), unique(d$cl))
  # same result when clustering is defined by the estimate object
  a3 <- estimate(estimate(g, id = d$cl), f, data = d, id = d$cl,
                 average = TRUE)
  expect_equivalent(vcov(a3), vcov(a2))

  # ids of 'data' default to index(x) for estimate objects of matching length
  a4 <- estimate(estimate(g, id = d$id), f, data = d, average = TRUE)
  expect_equivalent(vcov(a4), vcov(a1))
  expect_identical(index(a4), d$id)
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
  e1b <- summary(estimate(a3d, f=B), null=null)
  e1a <- summary(estimate(e1, f=diag(3)), null=null)
  expect_equal(unname(df), unname(e1a$compare$parameter))
  expect_equal(unname(df), unname(e1b$compare$parameter))
  expect_true(abs(pchisq(q, df=df, lower.tail=FALSE) - e1a$compare$p.value)<1e-16)
  expect_true(abs(pchisq(q, df=df, lower.tail=FALSE) - e1b$compare$p.value)<1e-16)
  expect_true(abs(q - e1a$compare$statistic) < 1e-9)
  expect_true(abs(q - e1b$compare$statistic) < 1e-9)
  # function spec.
  e2 <- estimate(a3d, function(p) c(p[1]-p[2], p[2]-p[3], p[1]-p[3]))
  e2b <- summary(estimate(a3d, function(p) c(p[1]-p[2], p[2]-p[3], p[1]-p[3])), null=null)
  e2a <- summary(estimate(e2, f=diag(3)), null=null)
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
  e1b <- summary(estimate(a3d, f=B), null=null)
  e1a <- summary(estimate(e1, f=diag(4)), null=null)
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
  e    <- summary(estimate(a1, f = B), null = null)
  p <- coef(a1)
  S <- vcov(a1)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
})

test_that("1D: Identity contrast, null = 1", {
  B    <- matrix(1, nrow = 1, ncol = 1)
  null <- 1
  e <- summary(estimate(a1, f = B), null = null)
  p <- coef(a1)
  S <- vcov(a1)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
})

test_that("1D: Scalar contrast (scaling), null = 0", {
  B    <- matrix(2, nrow = 1, ncol = 1)
  null <- 0
  e <- summary(estimate(a1, f = B), null = null)
  p <- coef(a1)
  S <- vcov(a1)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
})

test_that("2D: Identity contrast, null = c(0, 0)", {
  B    <- diag(2)
  null <- c(0, 0)
  e <- summary(estimate(a, f = B), null=null)
  p <- coef(a)
  S <- vcov(a)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
})

test_that("2D: Identity contrast, null equals true coefs", {
  B    <- diag(2)
  null <- c(1, 2)
  e <- summary(estimate(a, f = B), null = null)
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
  e <- summary(estimate(a, f = B), null = null)
  p <- coef(a)
  S <- vcov(a)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
})

test_that("2D: Difference contrast c(1,-1), null = 0", {
  B    <- matrix(c(1, -1), nrow = 1)
  null <- 0
  e <- summary(estimate(a, f = B), null = null)
  p <- coef(a)
  S <- vcov(a)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
})

test_that("2D: Difference contrast c(1,-1), null = -1", {
  B    <- matrix(c(1, -1), nrow = 1)
  null <- -1  # Testing H0: a1 - a2 = -1 (true, since 1 - 2 = -1)
  e <- summary(estimate(a, f = B), null = null)
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
  e <- summary(estimate(a, f = B), null = null)
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
  e <- summary(estimate(a, f = B), null = null)
  p <- coef(a)
  S <- vcov(a)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
})

test_that("3D: Identity contrast, null = true coefs", {
  B    <- diag(3)
  null <- c(1, 2, 3)
  e <- summary(estimate(a3d, f = B), null = null)
  p <- coef(a3d)
  S <- vcov(a3d)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
  expect_lt(e$compare$statistic, 1e-10)
})

test_that("3D: Identity contrast, null = c(0, 0, 0)", {
  B    <- diag(3)
  null <- c(0, 0, 0)
  e <- summary(estimate(a3d, f = B), null = null)
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
  e <- summary(estimate(a3d, f = B), null = null)
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
  e <- summary(estimate(a3d, f = B), null = null)
  p <- coef(a3d)
  S <- vcov(a3d)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
  expect_lt(e$compare$statistic, 1e-10)
})

test_that("4D: Identity contrast, null equals true coefs", {
  B    <- diag(4)
  null <- c(1, 2, 3, 4)
  e <- summary(estimate(a4d, f = B), null = null)
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
  e <- summary(estimate(a4d, f = B), null = null)
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
  e <- summary(estimate(a4d, f = B), null = null)
  p <- coef(a4d)
  S <- vcov(a4d)
  expected <- compute_wald(B, p, S, null)
  expect_equivalent(e$compare$statistic, expected)
  expect_lt(e$compare$statistic, 1e-10)
})

test_that("Degrees of freedom equals nrow(B)", {
  B    <- diag(2)
  null <- c(0, 0)
  e    <- summary(estimate(a, f = B), null = null)
  expect_identical(unname(e$compare$parameter), nrow(B))
  B    <- matrix(c(1, -1, 0, 0, 1, -1), nrow = 2, byrow = TRUE)
  null <- c(0, 0)
  e    <- summary(estimate(a3d, f = B), null = null)
  expect_equal(unname(e$compare$parameter), nrow(B))
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
    estimate(coef = c(a = 1), IC = ic_bad)
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
  e <- summary(estimate(a3d, "a1", "a1"), null = c(0, 1))
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

test_that("initialization of object without names", {
  p <- c(a=1, b=2, c=3, d=4)
  e0 <- estimate(coef=p, vcov=diag(4))
  e1 <- estimate(coef=1:4, vcov=diag(4))
  expect_equal(coef(e0), coef(e1), check.attributes=FALSE)
  expect_equal(vcov(e0), vcov(e1), check.attributes=FALSE)

  e <- estimate(coef=1:4, vcov=diag(4), f=list(1,2,3))
  expect_equal(coef(e), coef(e0)[c(1:3)], check.attributes=FALSE)
})
