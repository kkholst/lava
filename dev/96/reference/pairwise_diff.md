# Pairwise parameter differences

Construct contrasts for pairwise differences between parameters, or
apply these contrasts to an `estimate` object.

## Usage

``` r
pairwise_diff(n, idx = NULL, reverse = FALSE)
```

## Arguments

- n:

  An integer giving the total number of parameters, or an `estimate`
  object from which the number of parameters is derived.

- idx:

  An integer vector of parameter indices to compare. If `NULL`, all
  parameters are compared.

- reverse:

  Logical; if `TRUE`, reverse the direction of each difference by
  multiplying the contrast matrix by -1.

## Value

A numeric contrast matrix if `n` is an integer, or an `estimate` object
if `n` is an `estimate` object.

## Details

When `n` is an integer, the function returns a contrast matrix with `n`
columns and one row for each pairwise comparison. If `idx` is supplied,
only pairs of parameters selected by `idx` are compared.

When `n` is an `estimate` object, the number of parameters is obtained
from `length(coef(n))`. The contrast matrix is applied to the object,
returning an `estimate` object containing the pairwise differences and
their propagated uncertainty.

## Examples

``` r
pairwise_diff(4)
#>      [,1] [,2] [,3] [,4]
#> [1,]    1   -1    0    0
#> [2,]    1    0   -1    0
#> [3,]    1    0    0   -1
#> [4,]    0    1   -1    0
#> [5,]    0    1    0   -1
#> [6,]    0    0    1   -1
pairwise_diff(4, idx = c(1, 3, 4))
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0   -1    0
#> [2,]    1    0    0   -1
#> [3,]    0    0    1   -1
pairwise_diff(4, reverse = TRUE)
#>      [,1] [,2] [,3] [,4]
#> [1,]   -1    1    0    0
#> [2,]   -1    0    1    0
#> [3,]   -1    0    0    1
#> [4,]    0   -1    1    0
#> [5,]    0   -1    0    1
#> [6,]    0    0   -1    1

e <- estimate(coef = c(a = 1, b = 2, c = 4),
              vcov = diag(3))
pairwise_diff(e)
#>           Estimate Std.Err   2.5%   97.5% P-value
#> [a] - [b]       -1   1.414 -3.772  1.7718 0.47950
#> [a] - [c]       -3   1.414 -5.772 -0.2282 0.03389
#> [b] - [c]       -2   1.414 -4.772  0.7718 0.15730
```
