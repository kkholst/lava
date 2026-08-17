# Split data into folds

Split data into folds

## Usage

``` r
csplit(x, p = NULL, replace = FALSE, return.index = FALSE, k = 2, ...)
```

## Arguments

- x:

  Data or integer (size)

- p:

  Number of folds, or if a number between 0 and 1 is given two folds of
  size p and (1-p) will be returned

- replace:

  With or with-out replacement

- return.index:

  If TRUE index of folds are returned otherwise the actual data splits
  are returned (default)

- k:

  (Optional, only used when p=NULL) number of folds without shuffling

- ...:

  additional arguments to lower-level functions

## Author

Klaus K. Holst

## Examples

``` r
foldr(5,2,rep=2)
#> [[1]]
#> [[1]]$`1`
#> [1] 4 1 3
#> 
#> [[1]]$`2`
#> [1] 2 5
#> 
#> 
#> [[2]]
#> [[2]]$`1`
#> [1] 1 5 4
#> 
#> [[2]]$`2`
#> [1] 3 2
#> 
#> 
csplit(10,3)
#> $`1`
#> [1] 6 7 8 4
#> 
#> $`2`
#> [1]  5  9 10
#> 
#> $`3`
#> [1] 2 3 1
#> 
csplit(iris[1:10,]) ## Split in two sets 1:(n/2) and (n/2+1):n
#> $`1`
#>   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
#> 1          5.1         3.5          1.4         0.2  setosa
#> 2          4.9         3.0          1.4         0.2  setosa
#> 3          4.7         3.2          1.3         0.2  setosa
#> 4          4.6         3.1          1.5         0.2  setosa
#> 5          5.0         3.6          1.4         0.2  setosa
#> 
#> $`2`
#>    Sepal.Length Sepal.Width Petal.Length Petal.Width Species
#> 6           5.4         3.9          1.7         0.4  setosa
#> 7           4.6         3.4          1.4         0.3  setosa
#> 8           5.0         3.4          1.5         0.2  setosa
#> 9           4.4         2.9          1.4         0.2  setosa
#> 10          4.9         3.1          1.5         0.1  setosa
#> 
csplit(iris[1:10,],0.5)
#> [[1]]
#>    Sepal.Length Sepal.Width Petal.Length Petal.Width Species
#> 8           5.0         3.4          1.5         0.2  setosa
#> 7           4.6         3.4          1.4         0.3  setosa
#> 9           4.4         2.9          1.4         0.2  setosa
#> 2           4.9         3.0          1.4         0.2  setosa
#> 10          4.9         3.1          1.5         0.1  setosa
#> 
#> [[2]]
#>   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
#> 6          5.4         3.9          1.7         0.4  setosa
#> 3          4.7         3.2          1.3         0.2  setosa
#> 4          4.6         3.1          1.5         0.2  setosa
#> 1          5.1         3.5          1.4         0.2  setosa
#> 5          5.0         3.6          1.4         0.2  setosa
#> 
```
