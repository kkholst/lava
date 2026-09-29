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

  (optional) contrast matrix for the final Wald test. When supplied
  together with `null`, tests \\H_0: B\theta = b_0\\.

- ...:

  additional arguments passed to
  [contr](https://kkholst.github.io/lava/reference/contr.md).

- null:

  (optional) null hypothesis to test (default 0).

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

## Details

types of small-sample corrections:

- `"robust"` (default): no correction.

- `"df"`: applies \\n/(n-p)\\ correction (Mancl & DeRouen, 2001).

- `"mbn"`: Morel-Bokossa-Neerchal (2003) correction.

- `"hc3"`: leverage-adjusted HC3-type correction (blended with
  `var.adj`).

- `"hc4"`: Cribari-Neto (2004) leverage-adjusted correction.

The var.adj parameter controls the blending parameter for the HC3
leverage adjustment, by controls the weight between observation-level
empirical leverage and the average leverage \\p/n\\.

## See also

[`estimate.default()`](https://kkholst.github.io/lava/reference/estimate.default.md)
