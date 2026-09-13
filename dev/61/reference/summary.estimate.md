# Summary of estimate objects

Computes hypothesis tests, contrasts, and small-sample corrections for
an
[estimate](https://kkholst.github.io/lava/reference/estimate.default.md)
object. The arguments `null`, `contrast`, `type`, and `var.adj` were
previously available on
[`estimate.default()`](https://kkholst.github.io/lava/reference/estimate.default.md)
and have been moved here.

## Usage

``` r
# S3 method for class 'estimate'
summary(
  object,
  contrast,
  ...,
  null = 0,
  level = 0.95,
  type,
  var.adj = 0.25,
  df,
  transform = NULL,
  print = NULL
)
```

## Arguments

- object:

  an `estimate` object.

- contrast:

  (optional) contrast matrix for the final Wald test.

- ...:

  additional arguments passed to
  [contr](https://kkholst.github.io/lava/reference/contr.md).

- null:

  (optional) null hypothesis to test.

- level:

  level of confidence limits (default 0.95)

- type:

  type of small-sample correction. Requires the estimate to have been
  computed with `IC=TRUE` (the default).

- var.adj:

  variance adjustment parameter for small-sample correction. Requires
  the estimate to have been computed with `IC=TRUE` (the default).

- df:

  degrees of freedom for t-based inference (default: `NULL` for Gaussian
  approximation; when set, confidence intervals and p-values use the
  t-distribution with `df` degrees of freedom)

- transform:

  (optional) function applied to the point estimates and confidence
  interval bounds *after* inference is performed on the original scale.
  Useful for variance-stabilizing transformations, e.g., compute CIs on
  the `atanh` (Fisher z) scale and back-transform with `tanh`.

- print:

  (optional) custom print function for the resulting `summary.estimate`
  object

## See also

[`estimate.default()`](https://kkholst.github.io/lava/reference/estimate.default.md)
