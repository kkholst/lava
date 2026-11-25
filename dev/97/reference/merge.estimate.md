# Merge estimate objects

Merge estimate objects

## Usage

``` r
# S3 method for class 'estimate'
merge(
  x,
  y,
  ...,
  id,
  paired = FALSE,
  labels = NULL,
  keep = NULL,
  subset = NULL,
  regex = FALSE,
  sep = FALSE,
  drop.ic = FALSE,
  ignore.case = FALSE,
  sort = FALSE
)
```

## Arguments

- x:

  Object of class `estimate`

- y:

  Object of class `estimate`

- ...:

  Additional `estimate` objects or arguments

- id:

  Optional cluster variable

- paired:

  If TRUE a paired (matched) analysis is performed

- labels:

  Optional character vector of labels for the merged estimates

- keep:

  Optional character vector of parameter names to keep

- subset:

  Optional character vector of parameter names to subset

- regex:

  If TRUE, `keep` and `subset` are treated as regular expressions

- sep:

  Separator used for labeling

- drop.ic:

  If TRUE, drop the influence function from the result

- ignore.case:

  If TRUE, case is ignored in `keep`/`subset` matching

- sort:

  if true the returned influence function will be sorted according to
  the id variables (lexigraphically)

## Value

Object of class `estimate` (see
[estimate.default](https://kkholst.github.io/lava/reference/estimate.default.md)).

## See also

[`c.estimate()`](https://kkholst.github.io/lava/reference/c.estimate.md)
