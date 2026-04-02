# Generalized matrix inverse

Generalized matrix inverse

## Usage

``` r
Inverse(
  X,
  tol = lava.options()$itol,
  det = TRUE,
  names = !chol,
  chol = FALSE,
  symmetric = FALSE
)
```

## Arguments

- X:

  nxn matrix

- tol:

  tolerance for pseudo inverse

- det:

  logical, if true the determinant is returned

- names:

  preserve dimnames

- chol:

  use Cholesky decomposition for calculating inverse otherwise SVD

- symmetric:

  set to true if matrix is symmetric
