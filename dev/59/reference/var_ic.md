# Variance based on influence function

Computes the empirical variance of the influence function, providing a
sandwich-type variance estimator.

## Usage

``` r
var_ic(x, ...)
```

## Arguments

- x:

  Influence function matrix (observations x parameters), typically
  obtained from
  [IC](https://kkholst.github.io/lava/reference/IC.default.md).

- ...:

  Additional arguments (currently not used)

## Value

Covariance matrix (p x p).
