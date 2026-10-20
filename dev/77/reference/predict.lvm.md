# Prediction in structural equation models

Prediction in structural equation models

## Usage

``` r
# S3 method for class 'lvm'
predict(
  object,
  x = NULL,
  y = NULL,
  residual = FALSE,
  p,
  data,
  path = FALSE,
  quick = is.null(x) & !(residual | path),
  ...
)
```

## Arguments

- object:

  Model object

- x:

  optional list of (endogenous) variables to condition on

- y:

  optional subset of variables to predict

- residual:

  If true the residuals are predicted

- p:

  Parameter vector

- data:

  Data to use in prediction

- path:

  Path prediction

- quick:

  If TRUE the conditional mean and variance given covariates are
  returned (and all other calculations skipped)

- ...:

  Additional arguments to lower level function

## See also

predictlvm

## Examples

``` r
m <- lvm(list(c(y1,y2,y3)~u,u~x)); latent(m) <- ~u
d <- sim(m,100)
e <- estimate(m,d)

## Conditional mean (and variance as attribute) given covariates
r <- predict(e)
## Best linear unbiased predictor (BLUP)
r <- predict(e,vars(e))
##  Conditional mean of y3 giving covariates and y1,y2
r <- predict(e,y3~y1+y2)
##  Conditional mean  gives covariates and y1
r <- predict(e,~y1)
##  Predicted residuals (conditional on all observed variables)
r <- predict(e,vars(e),residual=TRUE)
```
