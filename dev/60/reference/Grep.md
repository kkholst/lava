# Finds elements in vector or column-names in data.frame/matrix

Pattern matching in a vector or column names of a data.frame or matrix.

## Usage

``` r
Grep(x, pattern, subset = TRUE, ignore.case = TRUE, ...)
```

## Arguments

- x:

  vector, matrix or data.frame.

- pattern:

  regular expression to search for

- subset:

  If TRUE returns subset of data.frame/matrix otherwise just the
  matching column names

- ignore.case:

  Default ignore case

- ...:

  Additional arguments to 'grep'

## Value

A data.frame with 2 columns with the indices in the first and the
matching names in the second.

## See also

[`grep`](https://rdrr.io/r/base/grep.html), and
[`agrep`](https://rdrr.io/r/base/agrep.html) for approximate string
matching,

## Author

Klaus K. Holst

## Examples

``` r
data(iris)
head(Grep(iris,"(len)|(sp)"))
#>   Sepal.Length Petal.Length Species
#> 1          5.1          1.4  setosa
#> 2          4.9          1.4  setosa
#> 3          4.7          1.3  setosa
#> 4          4.6          1.5  setosa
#> 5          5.0          1.4  setosa
#> 6          5.4          1.7  setosa
```
