# Subset parameter of estimate object

Methods for extracting or sorting parameter estimates

## Usage

``` r
# S3 method for class 'estimate'
subset(x, keep, ...)

# S3 method for class 'estimate'
rev(x)

# S3 method for class 'estimate'
tail(x, n = 1, ...)

# S3 method for class 'estimate'
head(x, n = 1, ...)

# S3 method for class 'estimate'
sort(x, decreasing = FALSE, by = c("labels", "values"), ...)

# S3 method for class 'estimate'
x[i, ...]

# S3 method for class 'estimate'
with(data, expr, ...)
```

## Arguments

- x:

  `estimate` object

- keep:

  subset (integer or character vector) of parameters to keep

- ...:

  arguments passed to [utils::tail](https://rdrr.io/r/utils/head.html),
  [utils::head](https://rdrr.io/r/utils/head.html),
  [base::sort](https://rdrr.io/r/base/sort.html),
  [base::subset](https://rdrr.io/r/base/subset.html)

- by:

  sorted by "values" or "labels" of parameters

- i:

  Indices specifying elements to extract

- data:

  data to use for constructing an environment. For the default `with`
  method this may be an environment, a list, a data frame, or an integer
  as in `sys.call`. For `within`, it can be a list or a data frame.

- expr:

  expression to evaluate; particularly for
  [`within()`](https://rdrr.io/r/base/with.html) often a “compound”
  expression, i.e., of the form

         {
           a <- somefun()
           b <- otherfun()
           .....
           rm(unused1, temp)
         }

## See also

[utils::tail](https://rdrr.io/r/utils/head.html)
[base::rev](https://rdrr.io/r/base/rev.html)
[utils::head](https://rdrr.io/r/utils/head.html)
[base::sort](https://rdrr.io/r/base/sort.html)
[base::subset](https://rdrr.io/r/base/subset.html)

## Examples

``` r
e <- estimate(coef=4:1, IC=matrix(0,nrow=2,ncol=4), labels=letters[1:4])
subset(e, c(1,4))
#>   Estimate Std.Err 2.5% 97.5% P-value
#> a        4       0    4     4       0
#> d        1       0    1     1       0
sort(e, by="values")
#>   Estimate Std.Err 2.5% 97.5% P-value
#> d        1       0    1     1       0
#> c        2       0    2     2       0
#> b        3       0    3     3       0
#> a        4       0    4     4       0
head(e, 2)
#>   Estimate Std.Err 2.5% 97.5% P-value
#> a        4       0    4     4       0
#> b        3       0    3     3       0
rev(e)
#>   Estimate Std.Err 2.5% 97.5% P-value
#> d        1       0    1     1       0
#> c        2       0    2     2       0
#> b        3       0    3     3       0
#> a        4       0    4     4       0
```
