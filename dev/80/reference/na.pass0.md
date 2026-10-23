# Handle Missing Values in Objects

Returns the object with missing data replaced by zeros. This is
sometimes useful for example when working with inverse probability
weighting of the complete-case data.

## Usage

``` r
na.pass0(object, na.action = na.omit, row.wise = FALSE, value = 0, ...)
```

## Arguments

- object:

  data.frame, vector

- na.action:

  function used to identify missing values

- row.wise:

  when true entire row is set to zero whenever a missing data is
  observed in that row

- value:

  value to replace NAs with (default 0)

- ...:

  additional arguments to lower level functions

## See also

[`na.pass()`](https://rdrr.io/r/stats/na.fail.html),
[`na.omit()`](https://rdrr.io/r/stats/na.fail.html),
[`na.fail()`](https://rdrr.io/r/stats/na.fail.html)

## Examples

``` r
d <- data.frame(y=c(1,1,NA,2,NA,2), r=c(1,1,0,1,1,1))
na.pass0(d)
#>   y r
#> 1 1 1
#> 2 1 1
#> 3 0 0
#> 4 2 1
#> 5 0 1
#> 6 2 1
glm(y ~ 1, weights=d$r, data=d, na.action=na.pass0)
#> 
#> Call:  glm(formula = y ~ 1, data = d, weights = d$r, na.action = na.pass0)
#> 
#> Coefficients:
#> (Intercept)  
#>         1.2  
#> 
#> Degrees of Freedom: 4 Total (i.e. Null);  4 Residual
#> Null Deviance:       2.8 
#> Residual Deviance: 2.8   AIC: Inf
```
