# Performs a rotation in the plane

Performs a rotation in the plane

## Usage

``` r
rotate2(x, theta = pi)
```

## Arguments

- x:

  Matrix to be rotated (2 times n)

- theta:

  Rotation in radians

## Value

Returns a matrix of the same dimension as `x`

## Author

Klaus K. Holst

## Examples

``` r
rotate2(cbind(c(1,2),c(2,1)))
#>      [,1] [,2]
#> [1,]   -1   -2
#> [2,]   -2   -1
```
