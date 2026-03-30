# Matching operator (x not in y) oposed to the `%in%`-operator (x in y)

Matching operator

## Usage

``` r
x %ni% y
```

## Arguments

- x:

  vector

- y:

  vector of same type as `x`

## Value

A logical vector.

## See also

[`match`](https://rdrr.io/r/base/match.html)

## Author

Klaus K. Holst

## Examples

``` r
1:10 %ni% c(1,5,10)
#>  [1] FALSE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE FALSE
```
