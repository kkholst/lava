# Estimate method for formulas

Fit a GLM model specified by a formula and apply
[estimate.default](https://kkholst.github.io/lava/reference/estimate.default.md)
for robust inference.

## Usage

``` r
# S3 method for class 'formula'
estimate(
  x,
  data,
  weights,
  family = stats::gaussian,
  ...,
  model = "glm",
  lvm = FALSE
)
```

## Arguments

- x:

  formula specifying the model

- data:

  data.frame

- weights:

  optional weights

- family:

  GLM family (default `gaussian`)

- ...:

  Additional arguments to
  [estimate.default](https://kkholst.github.io/lava/reference/estimate.default.md)

- model:

  character string specifying the model type (default `"glm"`)

- lvm:

  logical, if TRUE use `estimate.lvm` instead

## Value

Object of class `estimate` (see
[estimate.default](https://kkholst.github.io/lava/reference/estimate.default.md)).
