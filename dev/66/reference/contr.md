# Create contrast matrix

Create contrast matrix typically for use with 'estimate' (Wald tests).

## Usage

``` r
contr(p, n, diff = TRUE, ...)
```

## Arguments

- p:

  index of non-zero entries (see example)

- n:

  Total number of parameters (if omitted the max number in p will be
  used)

- diff:

  If FALSE all non-zero entries are +1, otherwise the second non-zero
  element in each row will be -1.

- ...:

  Additional arguments to lower level functions

## Examples

``` r
contr(2,n=5)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    1    0    0    0
contr(as.list(2:4),n=5)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    1    0    0    0
#> [2,]    0    0    1    0    0
#> [3,]    0    0    0    1    0
contr(list(1,2,4),n=5)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    0    0    0    0
#> [2,]    0    1    0    0    0
#> [3,]    0    0    0    1    0
contr(c(2,3,4),n=5)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    1   -1    0    0
#> [2,]    0    1    0   -1    0
contr(list(c(1,3),c(2,4)),n=5)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    0   -1    0    0
#> [2,]    0    1    0   -1    0
contr(list(c(1,3),c(2,4),5))
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    0   -1    0    0
#> [2,]    0    1    0   -1    0
#> [3,]    0    0    0    0    1

parsedesign(c("aa","b","c"),"?","?",diff=c(FALSE,TRUE))
#>      [,1] [,2] [,3]
#> [1,]    0    1    0
#> [2,]    0    0    1
#> [3,]    0    1   -1

## All pairs comparisons:
pdiff <- function(n) lava::contr(lapply(seq(n-1), function(x) seq(x, n)))
pdiff(4)
#>      [,1] [,2] [,3] [,4]
#> [1,]    1   -1    0    0
#> [2,]    1    0   -1    0
#> [3,]    1    0    0   -1
#> [4,]    0    1   -1    0
#> [5,]    0    1    0   -1
#> [6,]    0    0    1   -1
```
