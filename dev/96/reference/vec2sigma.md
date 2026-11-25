# Covariance matrix from upper-triangular elements

Conversion between covariance matrices and their upper-triangular
elements

## Usage

``` r
vec2sigma(x, groups = NULL, visits = NULL, simplify = FALSE)

sigma2vec(x)
```

## Arguments

- x:

  (numeric) for `vec2sigma`, the vector of upper-triangular elements of
  \\G\\ \\p\times p\\ matrices (length \\G p (p + 1) / 2\\). For
  `sigma2vec`, a matrix or a (named) list of matrices.

- groups:

  (character) group names of the covariance matrices. If `NULL`, the
  groups are derived from the names of `x`, which must be of the form
  `<group><index>` (e.g., `grpA1, ..., grpAk, grpB1, ..., grpBk`).

- visits:

  (character) optional row and column names of the matrices.

- simplify:

  (logical) if `TRUE` and there is only a single group, the matrix is
  returned instead of a list.

## Value

`vec2sigma`: list of matrices named by group (or a matrix if
`simplify = TRUE` and there is a single group). `sigma2vec`: numeric
vector with the upper-triangular elements of each matrix.

## Details

`vec2sigma` constructs symmetric (covariance) matrices from a vector of
upper-triangular elements (column-major, including the diagonal), e.g.,
the covariance parameters returned by `estimate(fit, sigma = TRUE)` for
an mmrm model. The vector may contain the elements of several matrices
of the same dimension (one per group), stacked after each other.
`sigma2vec` is the inverse operation.

## Examples

``` r
S <- matrix(c(2, 1, 1, 3), 2)
x <- sigma2vec(list(sigma = S))
x
#> sigma1 sigma2 sigma3 
#>      2      1      3 
vec2sigma(x)
#> $sigma
#>      [,1] [,2]
#> [1,]    2    1
#> [2,]    1    3
#> 
vec2sigma(x, groups = "sigma", visits = c("v1", "v2"), simplify = TRUE)
#>    v1 v2
#> v1  2  1
#> v2  1  3

## Two groups
vec2sigma(c(1, 0.5, 2, 3, 1, 4), groups = c("A", "B"))
#> $A
#>      [,1] [,2]
#> [1,]  1.0  0.5
#> [2,]  0.5  2.0
#> 
#> $B
#>      [,1] [,2]
#> [1,]    3    1
#> [2,]    1    4
#> 
```
