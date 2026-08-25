# Cross-validated two-stage estimator

Cross-validated two-stage estimator for non-linear SEM

## Usage

``` r
twostageCV(
  model1,
  model2,
  data,
  control1 = list(trace = 0),
  control2 = list(trace = 0),
  knots.boundary,
  nmix = 1:4,
  df = 1:9,
  fix = TRUE,
  std.err = TRUE,
  nfolds = 5,
  rep = 1,
  messages = 0,
  ...
)
```

## Arguments

- model1:

  model 1 (exposure measurement error model)

- model2:

  model 2

- data:

  data.frame

- control1:

  optimization parameters for model 1

- control2:

  optimization parameters for model 1

- knots.boundary:

  boundary points for natural cubic spline basis

- nmix:

  number of mixture components

- df:

  spline degrees of freedom

- fix:

  automatically fix parameters for identification (TRUE)

- std.err:

  calculation of standard errors (TRUE)

- nfolds:

  Number of folds (cross-validation)

- rep:

  Number of repeats of cross-validation

- messages:

  print information (\>0)

- ...:

  additional arguments to lower

## Examples

``` r
if (FALSE)  ## Reduce Ex.Timings
m1 <- lvm( x1+x2+x3 ~ u, latent= ~u)
m2 <- lvm( y ~ 1 )
m <- functional(merge(m1,m2), y ~ u, value=function(x) sin(x)+x)
#> Error: object 'm1' not found
distribution(m, ~u1) <- dist_uniform(-6, 6)
#> Error: object 'm' not found
d <- sim(m,n=500,seed=1)
#> Error: object 'm' not found
nonlinear(m2) <- y~u1
if (requireNamespace('mets', quietly=TRUE)) {
  set.seed(1)
  val <- twostageCV(m1, m2, data=d, std.err=FALSE, df=2:6, nmix=1:2,
                  nfolds=2)
  val
}
#> Error: object 'm1' not found
 # \dontrun{} ## Reduce Ex.Timings
```
