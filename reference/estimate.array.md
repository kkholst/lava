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

  target parameter ("mean", "variance", "quantile")

- probs:

  numeric vector of probabilities (for type="quantile")

- ...:

  Additional arguments to lower level functions (i.e.,
  stats::density.default when type="quantile")
