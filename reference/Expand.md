# Create a Data Frame from All Combinations of Factors

Create a Data Frame from All Combinations of Factors

## Usage

``` r
Expand(`_data`, ...)
```

## Arguments

- \_data:

  Data.frame

- ...:

  vectors, factors or a list containing these

## Details

Simple wrapper of the 'expand.grid' function. If x is a table then a
data frame is returned with one row pr individual observation.

## Author

Klaus K. Holst

## Examples

``` r
dd <- Expand(iris, Sepal.Length=2:8, Species=c("virginica","setosa"))
summary(dd)
#>   Sepal.Length        Species 
#>  Min.   :2.00   setosa    :7  
#>  1st Qu.:3.25   versicolor:0  
#>  Median :5.00   virginica :7  
#>  Mean   :5.00                 
#>  3rd Qu.:6.75                 
#>  Max.   :8.00                 

T <- with(warpbreaks, table(wool, tension))
Expand(T)
#>     wool tension
#> 1      A       L
#> 1.1    A       L
#> 1.2    A       L
#> 1.3    A       L
#> 1.4    A       L
#> 1.5    A       L
#> 1.6    A       L
#> 1.7    A       L
#> 1.8    A       L
#> 2      B       L
#> 2.1    B       L
#> 2.2    B       L
#> 2.3    B       L
#> 2.4    B       L
#> 2.5    B       L
#> 2.6    B       L
#> 2.7    B       L
#> 2.8    B       L
#> 3      A       M
#> 3.1    A       M
#> 3.2    A       M
#> 3.3    A       M
#> 3.4    A       M
#> 3.5    A       M
#> 3.6    A       M
#> 3.7    A       M
#> 3.8    A       M
#> 4      B       M
#> 4.1    B       M
#> 4.2    B       M
#> 4.3    B       M
#> 4.4    B       M
#> 4.5    B       M
#> 4.6    B       M
#> 4.7    B       M
#> 4.8    B       M
#> 5      A       H
#> 5.1    A       H
#> 5.2    A       H
#> 5.3    A       H
#> 5.4    A       H
#> 5.5    A       H
#> 5.6    A       H
#> 5.7    A       H
#> 5.8    A       H
#> 6      B       H
#> 6.1    B       H
#> 6.2    B       H
#> 6.3    B       H
#> 6.4    B       H
#> 6.5    B       H
#> 6.6    B       H
#> 6.7    B       H
#> 6.8    B       H
```
