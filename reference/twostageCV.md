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
 ## Reduce Ex.Timings##'
m1 <- lvm( x1+x2+x3 ~ u, latent= ~u)
m2 <- lvm( y ~ 1 )
m <- functional(merge(m1,m2), y ~ u, value=function(x) sin(x)+x)
distribution(m, ~u1) <- uniform.lvm(-6,6)
d <- sim(m,n=500,seed=1)
nonlinear(m2) <- y~u1
if (requireNamespace('mets', quietly=TRUE)) {
  set.seed(1)
  val <- twostageCV(m1, m2, data=d, std.err=FALSE, df=2:6, nmix=1:2,
                  nfolds=2)
  val
}
#> ────────────────────────────────────────────────────────────────────────────────
#> Selected mixture model: 1 component
#>       AIC1
#> 1 5130.210
#> 2 5132.707
#> ────────────────────────────────────────────────────────────────────────────────
#> Selected spline model degrees of freedom: 3
#> Knots: -2.674 -0.7956 1.082 2.96 
#> 
#>      RMSE(nfolds=2, rep=1)
#> df:1              5.353550
#> df:2              5.260141
#> df:3              4.851035
#> df:4              5.329716
#> df:5              6.220957
#> df:6              5.792509
#> ────────────────────────────────────────────────────────────────────────────────
#> 
#>                     Estimate Std. Error Z-value P-value
#> Regressions:                                           
#>    y~u1_1            1.38092                           
#>    y~u1_2            0.02123                           
#>    y~u1_3           -0.08440                           
#> Intercepts:                                            
#>    y                -0.33435                           
#> Residual Variances:                                    
#>    y                 1.61964                           
```
