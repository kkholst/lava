# expectation (g-computation)

sim_expectation <- function(n = 300, seed = 1) {
  set.seed(seed)
  dat <- data.frame(w1 = rnorm(n),
                    a = rbinom(n, 1, 0.5),
                    z = rbinom(n, 1, 0.5))
  dat$y <- rbinom(n, 1, plogis(-0.5 + dat$w1 + dat$a))
  dat$id <- paste0("a", seq_len(n))
  dat
}
dat <- sim_expectation()

test_that("expectation agrees with estimate(..., average=TRUE)", {
  g <- glm(y ~ w1 + a, data = dat, family = binomial)
  f <- function(p, data) plogis(p[1] + p["w1"] * data[, "w1"] + p["a"])
  m <- expectation(f, object=g)
  a <- suppressWarnings(estimate(g, f, average = TRUE))
  expect_equal(coef(m), coef(a))
  expect_equal(vcov(m), vcov(a))
  expect_equal(IC(m), IC(a), ignore_attr = TRUE)
  expect_equal(index(m), rownames(dat), ignore_attr = TRUE)

  # manual influence function
  X <- cbind(1, dat$w1, 1)
  q <- as.vector(plogis(X %*% coef(g)))
  ic <- q - mean(q) + IC(g) %*% colMeans(X * q * (1 - q))
  expect_equal(unname(coef(m)), mean(q), ignore_attr = TRUE)
  expect_true(mean(IC(m)-ic)<1e-9)

  # multiple parameters and labels
  f2 <- function(p, data) {
    list(p0 = plogis(p[1] + p["w1"] * data[, "w1"]),
         p1 = plogis(p[1] + p["w1"] * data[, "w1"] + p["a"]))
  }
  m2 <- expectation(f2, object=g)
  expect_identical(names(coef(m2)), c("p0", "p1"))
  expect_equal(coef(m2)[2], coef(m), ignore_attr = TRUE)
  expect_equal(vcov(m2)[2, 2], vcov(m), ignore_attr = TRUE)
  m3 <- expectation(f2, object=g, labels = c("y0", "y1"))
  expect_identical(names(coef(m3)), c("y0", "y1"))
})

test_that("expectation: ids of the g-computation and model part", {
  n <- nrow(dat)
  dat1 <- subset(dat, z == 1)
  g <- glm(y ~ w1 + a, data = dat1, family = binomial)
  e <- estimate(g, id = dat1$id)
  f <- function(p, data) plogis(p[1] + p["w1"] * data[, "w1"] + p["a"])

  X <- cbind(1, dat$w1, 1)
  q <- as.vector(plogis(X %*% coef(g)))
  D <- colMeans(X * q * (1 - q))
  ic2 <- rep(0, n)
  ic2[dat$z == 1] <- IC(g) %*% D * n / nrow(dat1)
  ic <- q - mean(q) + ic2

  m <- expectation(f, object=e, data = dat, id = dat$id)
  expect_equal(unname(coef(m)), mean(q))
  expect_true(mean(IC(m)-ic)<1e-9)
  expect_identical(index(m), dat$id)

  # id as column name of data
  expect_equal(vcov(expectation(object=e, f, data = dat, id = "id")), vcov(m))

  # ids from the names of the values of 'f'
  fn <- function(p, data) {
    structure(f(p, data), names = data$id)
  }
  mn <- expectation(object=e, fn, data = dat)
  expect_equal(vcov(mn), vcov(m))
  expect_identical(index(mn), dat$id)
  # 'id' takes precedence over the names of 'f'
  expect_message(
    expectation(object=e, fn, data = dat, id = paste0("x", seq_len(n))),
    "independence"
  )

  # default: rownames(data) (overlapping with rownames of the model frame)
  m0 <- expectation(object=g, f, data = dat)
  expect_equal(vcov(m0), vcov(m))

  # id = NULL removes the index
  m1 <- expectation(object=g, f, data = dat, id = NULL)
  expect_null(index(m1))
  expect_null(rownames(IC(m1)))
  expect_equal(vcov(m1), vcov(m), ignore_attr = TRUE)
})

test_that("expectation: 'f' returning values for some rows of data", {
  dat <- sim_expectation()
  rownames(dat) <- paste0("d", seq_len(nrow(dat)))
  dat1 <- subset(dat, z == 1)
  g <- glm(y ~ w1 + a, data = dat1, family = binomial)
  e <- estimate(g) # ids: rownames(dat1)
  # values only for the rows with a == 1 (named by rownames(data))
  f <- function(p, data) {
    X <- model.matrix(~ w1 + a, data = subset(data, a == 1))
    plogis(X %*% p)
  }
  f_all <- function(p, data) plogis(model.matrix(~ w1 + a, data) %*% p)
  # default ids: names of the values of 'f'
  m <- expectation(object=e, f, data = dat)
  expect_equal(index(m)[seq_len(sum(dat$a == 1))],
                   rownames(dat)[dat$a == 1])
  expect_equal(unname(coef(m)), mean(f_all(coef(g), subset(dat, a == 1))))
  # same estimate and variance as the conditional average over all rows of
  # 'data' (the IF of the latter also contains the rows with a == 0)
  ms1 <- expectation(object=e, f_all, data = dat, subset = dat$a == 1)
  ms2 <- expectation(object=e, f_all, data = dat, subset = "a")
  ms3 <- expectation(object=e, f_all, data = dat, subset = ~a==1)
  # variables of a formula (not in 'data') are looked up in the environment
  # of the formula and not in the calling environment
  gen <- function(thr) ~ a == thr
  ms4 <- local({
    thr <- 0
    expectation(object = e, f_all, data = dat, subset = gen(1))
  })
  # columns of 'data' mask variables in the environment of the formula
  ms5 <- local({
    a <- rep(0, nrow(dat))
    expectation(object = e, f_all, data = dat, subset = ~ a == 1)
  })
  expect_equal(coef(m), coef(ms1))
  expect_equal(vcov(m), vcov(ms1))
  expect_equal(coef(m), coef(ms2))
  expect_equal(vcov(m), vcov(ms2))
  expect_equal(coef(m), coef(ms3))
  expect_equal(vcov(m), vcov(ms3))
  expect_equal(coef(m), coef(ms4))
  expect_equal(vcov(m), vcov(ms4))
  expect_equal(coef(m), coef(ms5))
  expect_equal(vcov(m), vcov(ms5))
  expect_identical(index(ms1), rownames(dat))
  # model ids partly outside data
  id1 <- dat1$id
  id1[1:100] <- seq_len(100)
  m1 <- expectation(object=estimate(g, id = id1), f_all, data = dat, id = "id",
                 subset = "a")
  expect_equal(nrow(IC(m1)), nrow(dat) + 100)
  # 'id' and 'subset' must have one value per value of 'f'
  expect_error(expectation(object=e, f, data = dat, id = "id"), "does not agree")
  expect_error(expectation(object=e, f, data = dat, id = dat$id), "does not agree")
  expect_error(expectation(object=e, f, data = dat, subset = ~w1 > 0), "does not agree")
  # unnamed values for some of the rows cannot be linked
  expect_error(
    expectation(object=e, function(p, data) as.vector(f(p, data)), data = dat),
    "Unable to link"
  )
})

test_that("expectation: id matched to the values of 'f' (no data)", {
  dat <- sim_expectation()
  rownames(dat) <- paste0("d", seq_len(nrow(dat)))
  dat1 <- subset(dat, z == 1)
  g <- glm(y ~ w1 + a, data = dat1, family = binomial)
  id1 <- dat1$id
  id1[1:100] <- seq_len(100)
  e <- estimate(g, id = id1)
  # 'f' does not use 'data': one value per row of 'dat'
  f <- function(p, object) {
    plogis(p[1] + p["w1"] * dat$w1 + p["a"])
  }
  m <- expectation(object=e, f, id = dat$id)
  expect_equal(nrow(IC(m)), nrow(dat) + 100)
  expect_identical(index(m)[seq_len(nrow(dat))], dat$id)
  f_d <- function(p, data) plogis(p[1] + p["w1"] * data$w1 + p["a"])
  expect_equal(vcov(m), vcov(expectation(object=e, f_d, data = dat, id = "id")))
  # ids from the names of the values of 'f'
  fn <- function(p, object) structure(f(p, object), names = dat$id)
  expect_equal(vcov(expectation(object=e, fn)), vcov(m))
})

test_that("expectation: subset (conditional average)", {
  dat <- sim_expectation()
  n <- nrow(dat)
  dat1 <- subset(dat, z == 1)
  g <- glm(y ~ w1 + a, data = dat1, family = binomial)
  e <- estimate(g, id = dat1$id)
  f <- function(p, data) plogis(p[1] + p["w1"] * data[, "w1"] + p["a"])

  X <- cbind(1, dat$w1, 1)
  q <- as.vector(plogis(X %*% coef(g)))
  s <- dat$w1 > 0
  phat <- mean(s)
  mm <- mean(q * s)
  D_s <- colMeans(X * q * (1 - q) * s)
  ic2_s <- rep(0, n)
  ic2_s[dat$z == 1] <- IC(g) %*% D_s * n / nrow(dat1)
  icc <- (q * s - mm + ic2_s) / phat - mm / phat^2 * (s - phat)

  m <- expectation(object=e, f, data = dat, id = dat$id, subset = s)
  expect_equal(coef(m), mm / phat, ignore_attr = TRUE)
  expect_equal(IC(m), icc, ignore_attr = TRUE)
  # expression evaluated in data, and column name
  expect_equal(vcov(expectation(object=e, f, data = dat,
                                id = "id", subset = ~w1 > 0)),
               vcov(m))
  dat$s <- s
  expect_equal(vcov(expectation(object=e, f, data = dat,
                                id = "id", subset = "s")),
               vcov(m))
  # via estimate
})

test_that("expectation: arguments of 'f'", {
  g <- glm(y ~ w1 + a, data = dat, family = binomial)
  f <- function(p, data) plogis(p[1] + p["w1"] * data[, "w1"] + p["a"])
  m <- expectation(object=g, f)
  # 'object' argument (no data needed)
  fo <- function(p, object) {
    X <- model.matrix(object)
    X[, "a"] <- 1
    plogis(X %*% p)
  }
  expect_equal(vcov(expectation(object=g, fo)), vcov(m))
  # 'object' is the model also for estimate objects of the model
  fm <- function(p, data, object) {
    if (is.null(data)) data <- model.frame(object$fit)
    plogis(model.matrix(object$fit, data = transform(data, a = 1)) %*% p)
  }
  expect_equal(vcov(expectation(object=estimate(g), fm)), vcov(m))
  # dots
  fd <- function(p, ...) {
    data <- list(...)$data
    plogis(p[1] + p["w1"] * data[, "w1"] + p["a"])
  }
  expect_equal(vcov(expectation(object=g, fd)), vcov(m))
  # additional arguments
  fa <- function(p, data, a) plogis(p[1] + p["w1"] * data[, "w1"] + a * p["a"])
  expect_equal(vcov(expectation(object=g, fa, a = 1)), vcov(m))

  expect_error(expectation(object=g, function(x) x[1]), "must have an argument")
  expect_error(expectation(object=estimate(g, vcov = vcov(g)), f, data = dat),
               "Influence function")
})

test_that("estimate standardization", {
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
  est <- lava::expectation(object=g, function(p,data) {
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
  testthat::expect_equal(D1, D1a, ignore_attr = TRUE)
  testthat::expect_equal(D0, D0a, ignore_attr = TRUE)
  ic1 <- q1 - mean(q1) + apply(lava::IC(g), 1, \(x) sum(x*D1a))
  ic0 <- q0 - mean(q0) + apply(lava::IC(g), 1, \(x) sum(x*D0a))
  est2 <- c(y1=lava::estimate(coef=mean(q1), IC=ic1),
            y0=lava::estimate(coef=mean(q0), IC=ic0))
  testthat::expect_equal(coef(est), coef(est2), ignore_attr = TRUE)
  testthat::expect_equal(vcov(est), vcov(est2), ignore_attr = TRUE)

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
  a <- expectation(object=e, f, data = dat, id = dat$id)
  expect_equal(coef(a), mean(q), ignore_attr = TRUE)
  expect_equal(vcov(a), sum(ic^2) / n^2, ignore_attr = TRUE)
  expect_equal(IC(a), ic, ignore_attr = TRUE)
  expect_identical(index(a), dat$id)

  # plain glm: model identified by the rownames of its model frame, which do
  # not overlap with dat$id (independence)
  expect_message(
    a2 <- expectation(object=g, f, data = dat, id = dat$id),
    "independence"
  )
  expect_equal(coef(a2), coef(a), ignore_attr = TRUE)
  expect_equal(
    vcov(a2),
    var_ic(q - mean(q)) + var_ic(IC(g) %*% D), # independence
    ignore_attr = TRUE
  )
  expect_identical(index(a2), c(dat$id, rownames(dat1)))
  a3 <- expectation(object=g, f, data = dat) # id defaults to rownames
  expect_equal(vcov(a3), vcov(a), ignore_attr = TRUE)

  # row order of 'data' does not matter
  ord <- sample(n)
  a4 <- expectation(object=e, f, data = dat[ord, ], id = dat$id[ord])
  expect_equal(vcov(a4), vcov(a), ignore_attr = TRUE)
  expect_equal(IC(a4)[dat$id, ], IC(a)[dat$id, ], ignore_attr = TRUE)

  # conditional average (subset)
  # target: E[U(W_1, A = 1) | S = 1], for subset indicator S
  s <- dat$w1 > 0
  phat <- mean(s)
  m <- mean(q * s)
  D_s <- colMeans(X * q * (1 - q) * s)
  ic2_s <- rep(0, n)
  ic2_s[dat$z == 1] <- IC(g) %*% D_s * n / nrow(dat1)
  icc <- (q * s - m + ic2_s) / phat - m / phat^2 * (s - phat)
  ac <- expectation(object=e, f, data = dat, id = dat$id, subset = s)
  expect_equal(coef(ac), m / phat, ignore_attr = TRUE)
  expect_equal(vcov(ac), sum(icc^2) / n^2, ignore_attr = TRUE)
  expect_equal(IC(ac)[,1], icc, ignore_attr = TRUE)
  expect_equal(ac$id, dat$id, ignore_attr = TRUE)

  # multiple parameters
  # target: E[U(W_1, A = 1)], different U
  f2 <- function(p, data) {
    list(p0 = plogis(p[1] + p["w1"] * data[, "w1"]),
         p1 = plogis(p[1] + p["w1"] * data[, "w1"] + p["a"]))
  }
  a5 <- expectation(object=e, f2, data = dat, id = dat$id)
  expect_equal(coef(a5)[2], coef(a), ignore_attr = TRUE)
  expect_equal(vcov(a5)[2, 2], vcov(a), ignore_attr = TRUE)
  expect_equal(IC(a5)[,2], IC(a)[,1], ignore_attr = TRUE)

  # disjoint ids: independence between model and new data
  e_b <- estimate(g, id = paste0("b", seq_len(nrow(dat1))))
  expect_message(
    a6 <- expectation(object=e_b, f, data = dat, id = dat$id),
    "independence"
  )
  expect_equal(
    vcov(a6),
    sum((q - mean(q))^2) / n^2 + var_ic(IC(g) %*% D), ignore_attr = TRUE
  )
  expect_equal(nrow(IC(a6)), n + nrow(dat1))
  expect_identical(index(a6), c(dat$id, index(e_b)))

  # ids of estimate objects are used as is (not mapped via rownames of data):
  # default (rowname) ids of estimate(g) do not overlap with dat$id
  v_indep <- sum((q - mean(q))^2) / n^2 + var_ic(IC(g) %*% D)
  expect_message(
    a7 <- expectation(object=estimate(g), f, data = dat, id = dat$id),
    "independence"
  )
  expect_equal(nrow(IC(a7)), n + nrow(dat1))
  expect_equal(vcov(a7), v_indep, ignore_attr = TRUE)

  # estimate object without index: rownames of the IF are the ids
  ic_g <- IC(g)
  e0 <- estimate(coef = coef(g), IC = ic_g)
  expect_identical(index(e0), rownames(ic_g))
  expect_message(
    a8 <- expectation(object=e0, f, data = dat, id = dat$id),
    "independence"
  )
  expect_equal(nrow(IC(a8)), n + nrow(dat1))
  expect_equal(vcov(a8), v_indep, ignore_attr = TRUE)

  # estimate object without any ids: no link to data of different size
  # (no rownames in IC either)
  expect_error(
    expectation(object=estimate(g, id = NULL), f, data = dat, id = dat$id),
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
  a <- expectation(object=estimate(g, id = dm$id), f, data = dd, id = dd$id)
  expect_equal(nrow(IC(a)), N)
  expect_identical(index(a), pop$id)

  # manual: a_i = 1(i in data)(f_i - psi)/n_d + 1(i in model) D phi_i / n_m
  X <- cbind(1, dd$w, 1)
  q <- as.vector(plogis(X %*% coef(g)))
  D <- colMeans(X * q * (1 - q))
  ai <- rep(0, N)
  ai[1:600] <- (q - mean(q)) / 600
  ai[401:1000] <- ai[401:1000] + IC(g) %*% D / 600
  expect_equal(coef(a), mean(q), ignore_attr = TRUE)
  expect_equal(vcov(a), sum(ai^2), ignore_attr = TRUE)
  expect_equal(IC(a), N * ai, ignore_attr = TRUE)
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
  a0 <- expectation(object=g, f)
  # 'id' only refers to the g-computation part (the values of 'f'): the ids of
  # the model (rownames of the model frame) must be supplied separately
  expect_message(
    expectation(object=g, f, id = d$id),
    "independence"
  )
  # (estimate objects have no model frame, so 'data' is needed)
  a1 <- expectation(object=estimate(g, id = d$id), f, data = d, id = d$id)
  expect_equal(vcov(a1), vcov(a0), ignore_attr = TRUE)
  expect_equal(IC(a1), IC(a0), ignore_attr = TRUE)
  expect_identical(index(a1), d$id)

  X <- model.matrix(g)
  q <- as.vector(plogis(X %*% coef(g)))
  ic <- q - mean(q) + IC(g) %*% colMeans(X * q * (1 - q))
  icc <- rowsum(ic, d$cl, reorder = FALSE) * 50 / n
  a2 <- expectation(object=estimate(g, id = d$cl), f, data = d, id = d$cl)
  expect_equal(vcov(a2), var_ic(icc), ignore_attr = TRUE)
  expect_identical(index(a2), unique(d$cl))
  # same result with 'id' as column name
  a3 <- expectation(object=estimate(g, id = d$cl), f, data = d, id = "cl")
  expect_equal(vcov(a3), vcov(a2), ignore_attr = TRUE)

  # ids of 'data' default to rownames(data) (index(x) is not used)
  expect_message(
    a4 <- expectation(object=estimate(g, id = d$id), f, data = d),
    "independence"
  )
  expect_identical(index(a4), c(rownames(d), d$id))
  a5 <- expectation(object=estimate(g, id = d$id), f, data = d, id = "id")
  expect_equal(vcov(a5), vcov(a1), ignore_attr = TRUE)
  expect_identical(index(a5), d$id)
})

test_that("expectation: expressions and quoted arguments", {
  d <- data.frame(z=rnorm(10))
  e1 <- expression(z)
  e2 <- quote(z)
  m1 <- expectation(z, data=d)
  m2 <- expectation(e1, data=d)
  m3 <- expectation(e2, data=d)
  m4 <- with(d, expectation(z))
  m5 <- expectation(d$z)

  comp <- function(m1, m2) {
    expect_equal(coef(m1), coef(m2), ignore_attr=TRUE)
    expect_equal(IC(m1), IC(m2), ignore_attr=TRUE)
  }
  for (m in list(m2, m3, m4, m5)) comp(m1, m)
})
