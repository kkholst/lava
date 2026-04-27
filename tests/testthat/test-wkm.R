context("Weighted K-means (wkm)")

## Helper: best-permutation accuracy of cluster labels vs truth
.cluster_accuracy <- function(pred, truth) {
    tab <- table(pred, truth)
    K <- nrow(tab)
    perms_list <- function(k) {
        if (k == 1L) return(list(1L))
        sub <- perms_list(k - 1L)
        out <- list()
        for (s in sub) {
            for (i in seq_len(k)) {
                out[[length(out) + 1L]] <- append(s, k, after = i - 1L)
            }
        }
        out
    }
    perms <- perms_list(K)
    n_total <- sum(tab)
    best <- 0
    for (p in perms) {
        ## p is a permutation of columns (truth) to align with rows (pred)
        acc <- sum(diag(tab[, p, drop = FALSE])) / n_total
        if (acc > best) best <- acc
    }
    best
}

## ---------------------------------------------------------------------------
## Return structure
## ---------------------------------------------------------------------------
test_that("wkm returns list with cluster, center, ssw of expected types", {
    set.seed(1)
    x <- matrix(rnorm(60), ncol = 2)
    res <- wkm(x, mu = 2, n.start = 1)
    testthat::expect_type(res, "list")
    testthat::expect_named(res, c("cluster", "center", "ssw"),
                           ignore.order = TRUE)
    testthat::expect_length(res$cluster, nrow(x))
    testthat::expect_true(inherits(res$center, "by"))
    testthat::expect_equal(length(res$center), 2L)
    testthat::expect_true(all(res$ssw >= 0))
})

## ---------------------------------------------------------------------------
## Two well-separated clusters: correctness of clustering and centers
## ---------------------------------------------------------------------------
test_that("wkm separates two well-separated Gaussian blobs", {
    set.seed(2)
    n <- 100
    x <- rbind(matrix(rnorm(n * 2, mean = -3), ncol = 2),
               matrix(rnorm(n * 2, mean =  3), ncol = 2))
    truth <- c(rep(1L, n), rep(2L, n))
    res <- wkm(x, mu = 2, n.start = 5)
    testthat::expect_equal(length(unique(res$cluster)), 2L)
    ## centers (in some order) close to (-3,-3) and (3,3)
    cents <- do.call(rbind, lapply(res$center, as.numeric))
    ## sort by first coordinate
    cents <- cents[order(cents[, 1]), ]
    testthat::expect_lt(max(abs(cents[1, ] - c(-3, -3))), 0.5)
    testthat::expect_lt(max(abs(cents[2, ] - c(3, 3))), 0.5)
    ## cluster purity >= 95%
    acc <- .cluster_accuracy(res$cluster, truth)
    testthat::expect_gte(acc, 0.95)
})

## ---------------------------------------------------------------------------
## Custom initial centers (random.start = FALSE branch)
## ---------------------------------------------------------------------------
test_that("wkm accepts list of initial centers and forces n.start = 1", {
    set.seed(3)
    n <- 80
    x <- rbind(matrix(rnorm(n * 2, mean = -3), ncol = 2),
               matrix(rnorm(n * 2, mean =  3), ncol = 2))
    truth <- c(rep(1L, n), rep(2L, n))
    res <- wkm(x, mu = list(c(-3, -3), c(3, 3)), n.start = 7)
    testthat::expect_equal(length(res$center), 2L)
    acc <- .cluster_accuracy(res$cluster, truth)
    testthat::expect_gte(acc, 0.95)
})

## ---------------------------------------------------------------------------
## Univariate input
## ---------------------------------------------------------------------------
test_that("wkm handles univariate input via cbind", {
    set.seed(4)
    x <- c(rnorm(40, -2), rnorm(40, 2))
    res <- wkm(x, mu = 2, n.start = 3)
    testthat::expect_length(res$cluster, length(x))
    testthat::expect_equal(length(res$center), 2L)
})

## ---------------------------------------------------------------------------
## Formula interface
## ---------------------------------------------------------------------------
test_that("wkm formula interface works with data argument", {
    set.seed(5)
    res <- wkm(~ Sepal.Length + Sepal.Width, data = iris, mu = 3, n.start = 3)
    testthat::expect_length(res$cluster, nrow(iris))
    testthat::expect_lte(length(unique(res$cluster)), 3L)
})

## ---------------------------------------------------------------------------
## Weights effect: heavily weighted observations dominate within-cluster mean
## ---------------------------------------------------------------------------
test_that("wkm weights influence cluster centers", {
    testthat::skip_on_cran()
    set.seed(6)
    n <- 100
    ## Two well-separated blobs; within blob 2 add a heavily weighted shift
    blob1 <- matrix(rnorm(n * 2, mean = -3), ncol = 2)
    blob2 <- matrix(rnorm(n * 2, mean =  3), ncol = 2)
    x <- rbind(blob1, blob2)

    ## Uniform weights -> blob 2 center near (3, 3)
    res_u <- wkm(x, mu = list(c(-3, -3), c(3, 3)),
                 weights = rep(1, nrow(x)), n.start = 1)
    cu2 <- as.numeric(res_u$center[[2]])

    ## Heavily up-weight the *first half* of blob 2 (which is centered above
    ## the blob mean by construction of rnorm + sort by index).  The new
    ## center should be measurably different from the uniform-weight center.
    w <- rep(1, nrow(x))
    w[(n + 1):(n + n %/% 2)] <- 50
    res_w <- wkm(x, mu = list(c(-3, -3), c(3, 3)),
                 weights = w, n.start = 1)
    cw2 <- as.numeric(res_w$center[[2]])

    ## The weighted center should differ noticeably from the uniform center
    testthat::expect_gt(sqrt(sum((cw2 - cu2)^2)), 0.05)
    ## And it should equal the weighted mean of blob 2 observations assigned
    ## to cluster 2 (sanity check on lines 71-72).
    cl2_idx <- which(res_w$cluster == 2)
    expected <- colSums(x[cl2_idx, , drop = FALSE] * w[cl2_idx]) / sum(w[cl2_idx])
    testthat::expect_equal(cw2, as.numeric(expected), tolerance = 1e-8)
})

## ---------------------------------------------------------------------------
## iter.max = 1 still returns valid output
## ---------------------------------------------------------------------------
test_that("wkm runs with iter.max = 1", {
    set.seed(7)
    x <- matrix(rnorm(40), ncol = 2)
    res <- wkm(x, mu = 2, iter.max = 1, n.start = 1)
    testthat::expect_length(res$cluster, nrow(x))
    testthat::expect_equal(length(res$center), 2L)
})

## ---------------------------------------------------------------------------
## init argument: kmpp and fallback
## ---------------------------------------------------------------------------
test_that("wkm accepts init='kmpp' and falls back when init does not exist", {
    set.seed(8)
    x <- matrix(rnorm(60), ncol = 2)
    res1 <- wkm(x, mu = 2, init = "kmpp", n.start = 1)
    testthat::expect_length(res1$cluster, nrow(x))
    set.seed(8)
    res2 <- wkm(x, mu = 2, init = "this_init_does_not_exist", n.start = 1)
    testthat::expect_length(res2$cluster, nrow(x))
    testthat::expect_equal(length(res2$center), 2L)
})

## ---------------------------------------------------------------------------
## Reproducibility under set.seed
## ---------------------------------------------------------------------------
test_that("wkm is reproducible under set.seed", {
    x <- matrix(rnorm(80), ncol = 2)
    set.seed(42)
    a <- wkm(x, mu = 2, n.start = 3)
    set.seed(42)
    b <- wkm(x, mu = 2, n.start = 3)
    testthat::expect_equal(a$cluster, b$cluster)
    testthat::expect_equal(a$ssw, b$ssw)
})
