# Estimate parameters and influence function.

Estimate parameters for the sample mean, variance, and quantiles

## Usage

``` r
# S3 method for class 'array'
estimate(x, type = "mean", probs = 0.5, ...)
```

## Arguments

- x:

  numeric matrix

- type:

  target parameter ("mean", "variance", "quantile"). The type of
  quantile (see
  [`stats::quantile`](https://rdrr.io/r/stats/quantile.html)) can be
  chosen with "quantile1", "quantile2", ...

- probs:

  numeric vector of probabilities (for type="quantile")

- ...:

  Additional arguments to lower level functions (i.e.,
  stats::density.default when type="quantile")

## Value

Object of class `estimate` (see
[estimate.default](https://kkholst.github.io/lava/reference/estimate.default.md)).
