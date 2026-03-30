# Stack estimating equations

Stack estimating equations (two-stage estimator)

## Usage

``` r
# S3 method for class 'estimate'
stack(
  x,
  model2,
  D1u,
  inv.D2u,
  propensity,
  dpropensity,
  U,
  keep1 = FALSE,
  propensity.arg,
  estimate.arg,
  na.action = na.pass,
  ...
)
```

## Arguments

- x:

  Model 1

- model2:

  Model 2

- D1u:

  Derivative of score of model 2 w.r.t. parameter vector of model 1

- inv.D2u:

  Inverse of deri

- propensity:

  propensity score (vector or function)

- dpropensity:

  derivative of propensity score wrt parameters of model 1

- U:

  Optional score function (model 2) as function of all parameters

- keep1:

  If FALSE only parameters of model 2 is returned

- propensity.arg:

  Arguments to propensity function

- estimate.arg:

  Arguments to 'estimate'

- na.action:

  Method for dealing with missing data in propensity score

- ...:

  Additional arguments to lower level functions

## Examples

``` r
m <- lvm(z0~x)
Missing(m, z ~ z0) <- r~x
distribution(m,~x) <- binomial.lvm()
p <- c(r=-1,'r~x'=0.5,'z0~x'=2)
beta <- p[3]/2
d <- sim(m,500,p=p,seed=1)
m1 <- estimate(r~x,data=d,family=binomial)
d$w <- d$r/predict(m1,type="response")
m2 <- estimate(z~1, weights=w, data=d)
(e <- stack(m1,m2,propensity=TRUE))
#>             Estimate Std.Err   2.5% 97.5%   P-value
#> (Intercept)   0.9076 0.08836 0.7344 1.081 9.454e-25
```
