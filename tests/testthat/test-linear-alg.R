context("Linear algebra functions")


test_that("Matrix operations:", {
    ## vec operator
    testthat::expect_equivalent(vec(diag(3)),c(1,0,0,0,1,0,0,0,1))
    testthat::expect_true(nrow(vec(diag(3),matrix=TRUE))==9)

    ## commutaion matrix
    A <- matrix(1:16 ,ncol=4)
    K <- commutation(A)
    testthat::expect_equivalent(K%*%as.vector(A),vec(t(A),matrix=TRUE))
})

test_that("blockdiag:", {
    ## Block diagonal
    A <- diag(3)+1
    B <- blockdiag(A,A,A,pad=NA)
    testthat::expect_equivalent(dim(B),c(9,9))
    testthat::expect_true(sum(is.na(B))==81-27)

    B <- matrix(1, ncol=2, nrow=2)
    res <- blockdiag(A, B, B)
    expect_equivalent(res[1:3, 1:3], A)
    expect_equivalent(res[4:5, 4:5], B)
    expect_equivalent(res[6:7, 6:7], B)
    expect_true(sum(res==0) == (7*7 - 3*3 - 2*2 - 2*2))
})

test_that("wrapvev", {
    testthat::expect_equivalent(wrapvec(5,2),c(3,4,5,1,2))
    testthat::expect_equivalent(wrapvec(seq(1:5),-1),c(5,1,2,3,4))
})

test_that("revdiag", {
    A <- revdiag(1:3)
    testthat::expect_equivalent(A,matrix(c(0,0,1,0,2,0,3,0,0),3))
    testthat::expect_equivalent(1:3,revdiag(A))
    revdiag(A) <- 4
    testthat::expect_equivalent(rep(4,3),revdiag(A))
    diag(A) <- 0
    offdiag(A) <- 5
    testthat::expect_true(sum(offdiag(A))==6*5)
})

test_that("Inv, matrix inverse", {
    A <- matrix(0,3,3)
    offdiag(A,type=3) <- 1:6
    B <- crossprod(A)

    testthat::expect_equivalent(solve(A),Inverse(A))
    testthat::expect_equivalent(det(B),attr(Inverse(B,chol=TRUE),"det"))
})

test_that("rotate2", {
  ##rotate2: default rotation by pi flips the point
  x <- cbind(c(1, 0), c(0, 1))  # 2x2 matrix (2 rows, 2 cols)
  result <- rotate2(x, theta = pi)
  expected <- cbind(c(-1, 0), c(0, -1))
  expect_equal(result, expected, tolerance = 1e-10)

  ##rotate2: rotation by 0 returns the same matrix
  x <- cbind(c(1, 2), c(3, 4))
  result <- rotate2(x, theta = 0)
  expect_equal(result, x, tolerance = 1e-10)

  ##rotate2: rotation by pi/2 rotates 90 degrees
  x <- rbind(c(1, 0))  # row vector (1x2)
  result <- rotate2(x, theta = pi / 2)
  expected <- rbind(c(0, -1))
  expect_equal(result, expected, tolerance = 1e-10)

  ##rotate2: two rotations of pi/2 equal one rotation of pi
  x <- cbind(c(1, 2), c(2, 1))
  result_double <- rotate2(rotate2(x, theta = pi / 2), theta = pi / 2)
  result_pi     <- rotate2(x, theta = pi)
  expect_equal(result_double, result_pi, tolerance = 1e-10)
})

test_that("rot2D", {
  ## returns a 2x2 matrix
  R <- rot2D(pi / 4)
  expect_equal(dim(R), c(2, 2))
  ## rot2D: rotation matrix is orthogonal (R^T * R = I)
  R <- rot2D(pi / 3)
  expect_equal(t(R) %*% R, diag(2), tolerance = 1e-10)
  ## rot2D: determinant of rotation matrix is 1
  R <- rot2D(pi / 6)
  expect_equal(det(R), 1, tolerance = 1e-10)
  ## rot2D: rotation by 0 returns identity matrix
  R <- rot2D(0)
  expect_equal(R, diag(2), tolerance = 1e-10)
})

test_that("rot3D", {
  ## returns a 3x3 matrix
  R <- rot3D(x = pi / 4, y = pi / 6, z = pi / 3)
  expect_equal(dim(R), c(3, 3))
  ## rot3D: no rotation returns identity matrix
  R <- rot3D()
  expect_equal(R, diag(3), tolerance = 1e-10)
  ## rot3D: rotation matrix is orthogonal (R^T * R = I)
  R <- rot3D(x = pi / 4, y = pi / 3, z = pi / 6)
  expect_equal(t(R) %*% R, diag(3), tolerance = 1e-10)
  ## rot3D: determinant of rotation matrix is 1
  R <- rot3D(x = pi / 5, y = pi / 7, z = pi / 9)
  expect_equal(det(R), 1, tolerance = 1e-10)
  ## rot3D: rotation around z-axis only leaves z-component unchanged
  R <- rot3D(z = pi / 2)
  v <- c(1, 0, 0)
  result <- R %*% v
  # x-axis rotated 90° around z should give (0, 1, 0)
  expect_equal(as.numeric(result), c(0, 1, 0), tolerance = 1e-10)
  ## rot3D: rotation around x-axis only leaves x-component unchanged
  R <- rot3D(x = pi / 2)
  v <- c(0, 1, 0)
  result <- R %*% v
  # y-axis rotated 90° around x should give (0, 0, 1)
  expect_equal(as.numeric(result), c(0, 0, 1), tolerance = 1e-10)
})
