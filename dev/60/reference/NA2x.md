# Convert to/from NA

Convert vector to/from NA

## Usage

``` r
NA2x(s, x = 0)
```

## Arguments

- s:

  The input vector (of arbitrary class)

- x:

  The elements to transform into `NA` resp. what to transform `NA` into.

## Value

A vector with same dimension and class as `s`.

## Author

Klaus K. Holst

## Examples

``` r
##' 
x2NA(1:10, 1:5)
#>  [1] NA NA NA NA NA  6  7  8  9 10
NA2x(x2NA(c(1:10),5),5)##'
#>  [1]  1  2  3  4  5  6  7  8  9 10
```
