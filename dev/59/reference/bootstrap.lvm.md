# Calculate bootstrap estimates of a lvm object

Draws non-parametric bootstrap samples

## Usage

``` r
# S3 method for class 'lvm'
bootstrap(x,R=100,data,fun=NULL,control=list(),
                          p, parametric=FALSE, bollenstine=FALSE,
                          constraints=TRUE,sd=FALSE, mc.cores,
                          future.args=list(future.seed=TRUE),
                          ...)

# S3 method for class 'lvmfit'
bootstrap(x,R=100,data=model.frame(x),
                             control=list(start=coef(x)),
                             p=coef(x), parametric=FALSE, bollenstine=FALSE,
                             estimator=x$estimator,weights=Weights(x),...)
```

## Arguments

- x:

  `lvm`-object.

- R:

  Number of bootstrap samples

- data:

  The data to resample from

- fun:

  Optional function of the (bootstrapped) model-fit defining the
  statistic of interest

- control:

  Options to the optimization routine

- p:

  Parameter vector of the null model for the parametric bootstrap

- parametric:

  If TRUE a parametric bootstrap is calculated. If FALSE a
  non-parametric (row-sampling) bootstrap is computed.

- bollenstine:

  Bollen-Stine transformation (non-parametric bootstrap) for bootstrap
  hypothesis testing.

- constraints:

  Logical indicating whether non-linear parameter constraints should be
  included in the bootstrap procedure

- sd:

  Logical indicating whether standard error estimates should be included
  in the bootstrap procedure

- mc.cores:

  Optional number of cores for parallel computing. If omitted
  future.apply will be used (see future::plan)

- future.args:

  arguments to future.apply::future_lapply

- ...:

  Additional arguments, e.g. choice of estimator.

- estimator:

  String definining estimator, e.g. 'gaussian' (see `estimator`)

- weights:

  Optional weights matrix used by `estimator`

## Value

A `bootstrap.lvm` object.

## See also

[`confint.lvmfit`](https://kkholst.github.io/lava/reference/confint.lvmfit.md)

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm(y~x)
d <- sim(m,100)
e <- estimate(lvm(y~x), data=d)
 ## Reduce Ex.Timings
B <- bootstrap(e,R=50,mc.cores=1)
B
#> Non-parametric bootstrap statistics (R=50):
#> 
#>      Estimate     Bias         Std.Err      2.5 %        97.5 %      
#> y    -0.012523283 -0.002513924  0.083715845 -0.175568919  0.123923908
#> y~x   1.007018982  0.000333004  0.073128648  0.888965542  1.154184192
#> y~~y  0.939676415 -0.014462789  0.110341543  0.746729951  1.146074112
#> 
```
