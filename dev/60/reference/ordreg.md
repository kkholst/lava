# Univariate cumulative link regression models

Ordinal regression models

## Usage

``` r
ordreg(
  formula,
  data = parent.frame(),
  offset,
  family = stats::binomial("probit"),
  start,
  fast = FALSE,
  ...
)
```

## Arguments

- formula:

  formula

- data:

  data.frame

- offset:

  offset

- family:

  family (default proportional odds)

- start:

  optional starting values

- fast:

  If TRUE standard errors etc. will not be calculated

- ...:

  Additional arguments to lower level functions

## Details

Let \\Y\in\\1,...,J\\\\ be the ordinal outcome and \$X\$ a vector of
covariates. The cumulative link model is given by \$\$ P(Y\leq j\|X=x) =
g(a_j - b^t x), j=1,...,J-1.\$\$ The default link function is the Probit
function, i.e. where \$g\$ is equal to the standard normal cumulative
distribution function. The proportional odds model is obtained with
`family=binomial(logit)`.

Note, the intercept parameters are parametrized such that they are
monotone increasing \\a_1 \< \cdots \< a\_{J-1}\\. To get the parameter
estimates of the actual \$a_j\$'s use the `summary` method.

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm(y~x)
ordinal(m,K=3) <- ~y
d <- sim(m,100)
e <- ordreg(y~x,d)
summary(e)
#> AIC:  145.3235 
#> 
#>     Estimate Std.Err    2.5%   97.5%   P-value
#> 0|1  -1.0227  0.1491 -1.3150 -0.7305 6.892e-12
#> 1|2  -0.5031  0.1441 -0.7855 -0.2207 4.798e-04
#> x     0.9674  0.1675  0.6391  1.2957 7.674e-09
#> 
```
