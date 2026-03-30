# Create random missing data

Generates missing entries in data.frame/matrix

## Usage

``` r
makemissing(
  data,
  p = 0.2,
  cols = seq_len(ncol(data)),
  rowwise = FALSE,
  nafun = function(x) x,
  seed = NULL
)
```

## Arguments

- data:

  data.frame

- p:

  Fraction of missing data in each column

- cols:

  Which columns (name or index) to alter

- rowwise:

  Should missing occur row-wise (either none or all selected columns are
  missing)

- nafun:

  (Optional) function to be applied on data.frame before return (e.g.
  `na.omit` to return complete-cases only)

- seed:

  Random seed

## Value

data.frame

## Author

Klaus K. Holst
