# Extract influence function from model object

Extract i.i.d. decomposition (influence function) from model object

## Usage

``` r
# Default S3 method
IC(x, bread, id = NULL, ...)
```

## Arguments

- x:

  model object

- bread:

  (optional) Inverse of derivative of mean score function

- id:

  (optional) id/cluster variable

- ...:

  additional arguments

## Value

Matrix with rows corresponding to observations and columns to
parameters, representing the estimated influence function. Attributes
include `bread` (derivative of the estimating equation).

## See also

[var_ic](https://kkholst.github.io/lava/reference/var_ic.md),
[iid](https://kkholst.github.io/lava/reference/iid.md)

## Examples

``` r
m <- lvm(y~x+z)
distribution(m, ~y+z) <- binomial.lvm("logit")
d <- sim(m,1e3)
g <- glm(y~x+z,data=d,family=binomial)
var_ic(IC(g))
#>               (Intercept)             x           z
#> (Intercept)  1.005602e-02 -4.173875e-05 -0.01007301
#> x           -4.173875e-05  6.989745e-03  0.00149390
#> z           -1.007301e-02  1.493900e-03  0.02068351
```
