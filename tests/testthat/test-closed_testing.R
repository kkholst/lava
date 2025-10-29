context("Multiple testing")

test_that("closed_testing", {
  m <- lvm()
  regression(m, c(y1, y2, y3, y4) ~ x) <- c(0, 0.25, 0, 0.25)
  regression(m, to=endogenous(m), from="u") <- 1
  variance(m, endogenous(m)) <- 1
  set.seed(1)
  d <- sim(m, 200)
  l1 <- lm(y1~x, d)
  l2 <- lm(y2~x, d)
  l3 <- lm(y3~x, d)
  l4 <- lm(y4~x, d)
  e <- merge(l1, l2, l3, l4, subset = 2)

  adj <- closed_testing(e)

  H1 <- list(
    h1234 = lava::compare(e, contrast = diag(4)),
    h123 = lava::compare(e, contrast = diag(4)[c(1, 2, 3), , drop = FALSE]),
    h124 = lava::compare(e, contrast = diag(4)[c(1, 2, 4), , drop = FALSE]),
    h134 = lava::compare(e, contrast = diag(4)[c(1, 3, 4), , drop = FALSE]),
    h12 = lava::compare(e, contrast = diag(4)[c(1, 2), , drop = FALSE]),
    h13 = lava::compare(e, contrast = diag(4)[c(1, 3), , drop = FALSE]),
    h14 = lava::compare(e, contrast = diag(4)[c(1, 4), , drop = FALSE]),
    h1 = lava::compare(e, contrast = diag(4)[c(1), , drop = FALSE])
  )

  expect_true(adj$raw.pval[[4]] == H1$h1234$p.value)

  expect_true(adj$p.value[1] ==
              max(lapply(H1, function(x) x$p.value) |> unlist()))

  H2 <- list(
    h1234 = lava::compare(e, contrast = diag(4)),
    h234 = lava::compare(e, contrast = diag(4)[c(2, 3, 4), , drop = FALSE]),
    h123 = lava::compare(e, contrast = diag(4)[c(1, 2, 3), , drop = FALSE]),
    h124 = lava::compare(e, contrast = diag(4)[c(1, 2, 4), , drop = FALSE]),
    h12 = lava::compare(e, contrast = diag(4)[c(1, 2), , drop = FALSE]),
    h23 = lava::compare(e, contrast = diag(4)[c(2, 3), , drop = FALSE]),
    h24 = lava::compare(e, contrast = diag(4)[c(2, 4), , drop = FALSE]),
    h2 = lava::compare(e, contrast = diag(4)[c(2), , drop = FALSE])
  )

  expect_true(adj$raw.pval[[4]] == H1$h1234$p.value)

  expect_true(adj$p.value[2] ==
              max(lapply(H2, function(x) x$p.value) |> unlist()))

})
