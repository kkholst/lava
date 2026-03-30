# Combine matrices to block diagonal structure

Combine matrices to block diagonal structure

## Usage

``` r
blockdiag(x, ..., pad = 0)
```

## Arguments

- x:

  Matrix

- ...:

  Additional matrices

- pad:

  Vyalue outside block-diagonal

## Author

Klaus K. Holst

## Examples

``` r
A <- diag(3)+1
blockdiag(A,A,A,pad=NA)
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    2    1    1   NA   NA   NA   NA   NA   NA
#>  [2,]    1    2    1   NA   NA   NA   NA   NA   NA
#>  [3,]    1    1    2   NA   NA   NA   NA   NA   NA
#>  [4,]   NA   NA   NA    2    1    1   NA   NA   NA
#>  [5,]   NA   NA   NA    1    2    1   NA   NA   NA
#>  [6,]   NA   NA   NA    1    1    2   NA   NA   NA
#>  [7,]   NA   NA   NA   NA   NA   NA    2    1    1
#>  [8,]   NA   NA   NA   NA   NA   NA    1    2    1
#>  [9,]   NA   NA   NA   NA   NA   NA    1    1    2
```
