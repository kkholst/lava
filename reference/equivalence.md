# Identify candidates of equivalent models

Identifies candidates of equivalent models

## Usage

``` r
equivalence(x, rel, tol = 0.001, k = 1, omitrel = TRUE, ...)
```

## Arguments

- x:

  `lvmfit`-object

- rel:

  Formula or character-vector specifying two variables to omit from the
  model and subsequently search for possible equivalent models

- tol:

  Define two models as empirical equivalent if the absolute difference
  in score test is less than `tol`

- k:

  Number of parameters to test simultaneously. For `equivalence` the
  number of additional associations to be added instead of `rel`.

- omitrel:

  if `k` greater than 1, this boolean defines wether to omit candidates
  containing `rel` from the output

- ...:

  Additional arguments to be passed to the lower-level functions

## See also

[`compare`](http://kkholst.github.io/lava/reference/compare.md),
[`modelsearch`](http://kkholst.github.io/lava/reference/modelsearch.md)

## Author

Klaus K. Holst
