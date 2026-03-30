# Trace operator

Calculates the trace of a square matrix.

## Usage

``` r
tr(x, ...)
```

## Arguments

- x:

  Square numeric matrix

- ...:

  Additional arguments to lower level functions

## Value

`numeric`

## See also

[`crossprod`](https://rdrr.io/r/base/crossprod.html),
[`tcrossprod`](https://rdrr.io/r/base/crossprod.html)

## Author

Klaus K. Holst

## Examples

``` r
tr(diag(1:5))
#> [1] 15
```
