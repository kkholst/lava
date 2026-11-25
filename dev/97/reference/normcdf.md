# Normal distribution

Density (`normpdf`), distribution function (`normcdf`), and quantile
function (`norminv`) of the normal distribution for objects of class
`estimate` (see
[estimate.default](https://kkholst.github.io/lava/reference/estimate.default.md))
in which case the delta method is applied to obtain the influence
function (and hence robust standard errors) of the transformed
parameter.

## Usage

``` r
normcdf(q, mean = 0, sd = 1, ...)

# Default S3 method
normcdf(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, ...)

normpdf(x, mean = 0, sd = 1, ...)

# Default S3 method
normpdf(x, mean = 0, sd = 1, log = FALSE, ...)

norminv(p, mean = 0, sd = 1, ...)

# Default S3 method
norminv(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, ...)
```

## Arguments

- mean:

  mean of the normal distribution (numeric or `estimate` object)

- sd:

  standard deviation of the normal distribution (numeric or `estimate`
  object)

- ...:

  additional arguments to lower level functions (e.g., `labels`)

- lower.tail:

  if `TRUE` (default) probabilities are \\P(X\leq x)\\, otherwise \\P(X
  \> x)\\

- x, q:

  vector of quantiles, or an `estimate` object

- log, log.p:

  if `TRUE` probabilities/densities are given as \\\log(p)\\

- p:

  vector of probabilities, or an `estimate` object

## Value

A numeric vector (default method) or an object of class `estimate`.

## Details

Any of the arguments `x`/`q`/`p`, `mean` and `sd` may be an `estimate`
object. When more than one of them is an `estimate` object the objects
are first merged (see
[merge.estimate](https://kkholst.github.io/lava/reference/merge.estimate.md))
and the joint delta method is applied, correctly propagating the
dependence between them.

## See also

[estimate.default](https://kkholst.github.io/lava/reference/estimate.default.md)
[stats::dnorm](https://rdrr.io/r/stats/Normal.html)

## Author

Klaus Kähler Holst

## Examples

``` r
e <- estimate(coef = c(a = 1, b = 2), vcov = diag(2) * 0.1)
normcdf(e)
#>   Estimate Std.Err   2.5%  97.5%   P-value
#> a   0.8413 0.07652 0.6914 0.9913 4.021e-28
#> b   0.9772 0.01707 0.9438 1.0107 0.000e+00
normpdf(e)
#>   Estimate Std.Err     2.5%  97.5%  P-value
#> a  0.24197 0.07652  0.09200 0.3919 0.001565
#> b  0.05399 0.03415 -0.01294 0.1209 0.113846
norminv(normcdf(e)) # identity
#>   Estimate Std.Err   2.5% 97.5%   P-value
#> a        1  0.3162 0.3802  1.62 1.565e-03
#> b        2  0.3162 1.3802  2.62 2.540e-10

## mean and sd may themselves be estimates
mu <- estimate(coef = c(mu = 0.5), vcov = 0.01)
normcdf(e, mean = mu, sd = 2)
#>   Estimate Std.Err 2.5% 97.5% P-value
#> a   0.5987                           
#> b   0.7734                           
```
