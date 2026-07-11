# Estimate method for GLM objects

Applies
[estimate.default](https://kkholst.github.io/lava/reference/estimate.default.md)
to obtain robust inference for a fitted `glm` object.

## Usage

``` r
# S3 method for class 'glm'
estimate(x, ...)
```

## Arguments

- x:

  Fitted `glm` object

- ...:

  Additional arguments to
  [estimate.default](https://kkholst.github.io/lava/reference/estimate.default.md)

## Value

Object of class `estimate` (see
[estimate.default](https://kkholst.github.io/lava/reference/estimate.default.md)).
