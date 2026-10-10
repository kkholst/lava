# Estimation of parameters in a Latent Variable Model (lvm)

Estimate parameters. MLE, IV or user-defined estimator.

## Usage

``` r
# S3 method for class 'lvm'
estimate(
  x,
  data = parent.frame(),
  estimator = NULL,
  control = list(),
  missing = FALSE,
  weights,
  weightsname,
  data2,
  id,
  fix,
  index = !quick,
  graph = FALSE,
  messages = lava.options()$messages,
  quick = FALSE,
  method,
  param,
  cluster,
  p,
  ...
)
```

## Arguments

- x:

  `lvm`-object

- data:

  `data.frame`

- estimator:

  String defining the estimator (see details below)

- control:

  control/optimization parameters (see details below)

- missing:

  Logical variable indiciating how to treat missing data. Setting to
  FALSE leads to complete case analysis. In the other case likelihood
  based inference is obtained by integrating out the missing data under
  assumption the assumption that data is missing at random (MAR).

- weights:

  Optional weights to used by the chosen estimator.

- weightsname:

  Weights names (variable names of the model) in case `weights` was
  given as a vector of column names of `data`

- data2:

  Optional additional dataset used by the chosen estimator.

- id:

  Vector (or name of column in `data`) that identifies correlated groups
  of observations in the data leading to variance estimates based on a
  sandwich estimator

- fix:

  Logical variable indicating whether parameter restriction
  automatically should be imposed (e.g. intercepts of latent variables
  set to 0 and at least one regression parameter of each measurement
  model fixed to ensure identifiability.)

- index:

  For internal use only

- graph:

  For internal use only

- messages:

  Control how much information should be printed during estimation (0:
  none)

- quick:

  If TRUE the parameter estimates are calculated but all additional
  information such as standard errors are skipped

- method:

  Optimization method

- param:

  set parametrization (see
  [`help(lava.options)`](https://kkholst.github.io/lava/reference/lava.options.md))

- cluster:

  Obsolete. Alias for 'id'.

- p:

  Evaluate model in parameter 'p' (no optimization)

- ...:

  Additional arguments to be passed to lower-level functions

## Value

A `lvmfit`-object.

## Details

A list of parameters controlling the estimation and optimization
procedures is parsed via the `control` argument. By default Maximum
Likelihood is used assuming multivariate normal distributed measurement
errors. A list with one or more of the following elements is expected:

- start::

  Starting value. The order of the parameters can be shown by calling
  `coef` (with `mean=TRUE`) on the `lvm`-object or with
  `plot(..., labels=TRUE)`. Note that this requires a check that it is
  actual the model being estimated, as `estimate` might add additional
  restriction to the model, e.g. through the `fix` and `exo.fix`
  arguments. The `lvm`-object of a fitted model can be extracted with
  the `Model`-function.

- starterfun::

  Starter-function with syntax `function(lvm, S, mu)`. Three builtin
  functions are available: `startvalues`, `startvalues0`,
  `startvalues1`, ...

- estimator::

  String defining which estimator to use (Defaults to “`gaussian`”)

- meanstructure:

  Logical variable indicating whether to fit model with meanstructure.

- method::

  String pointing to alternative optimizer (e.g. `optim` to use
  simulated annealing).

- control::

  Parameters passed to the optimizer (default
  [`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html)).

- tol::

  Tolerance of optimization constraints on lower limit of variance
  parameters.

## See also

estimate.default score, information

## Author

Klaus K. Holst

## Examples

``` r
dd <- read.table(header=TRUE,
text="x1 x2 x3
 0.0 -0.5 -2.5
-0.5 -2.0  0.0
 1.0  1.5  1.0
 0.0  0.5  0.0
-2.5 -1.5 -1.0")
e <- estimate(lvm(c(x1,x2,x3)~u),dd)

## Simulation example
m <- lvm(list(y~v1+v2+v3+v4,c(v1,v2,v3,v4)~x))
covariance(m) <- v1~v2+v3+v4
dd <- sim(m,10000) ## Simulate 10000 observations from model
e <- estimate(m, dd) ## Estimate parameters
e
#>                      Estimate Std. Error   Z-value  P-value
#> Regressions:                                               
#>    y~v1               1.00369    0.01783  56.30059   <1e-12
#>    y~v2               0.99186    0.01083  91.61548   <1e-12
#>    y~v3               1.00216    0.01101  91.01110   <1e-12
#>    y~v4               0.99616    0.01098  90.69100   <1e-12
#>     v1~x              0.99178    0.01011  98.07103   <1e-12
#>    v2~x               1.00144    0.01013  98.81021   <1e-12
#>     v3~x              0.98784    0.01025  96.33336   <1e-12
#>    v4~x               1.00687    0.00988 101.94449   <1e-12
#> Intercepts:                                                
#>    y                  0.00290    0.01000   0.29028   0.7716
#>    v1                -0.00637    0.01002  -0.63591   0.5248
#>    v2                 0.00609    0.01004   0.60632   0.5443
#>    v3                -0.00834    0.01016  -0.82076   0.4118
#>    v4                -0.01165    0.00979  -1.19001    0.234
#> Residual Variances:                                        
#>    y                  0.99981    0.01414  70.71068         
#>    v1                 1.00448    0.01122  89.53613         
#>    v1~~v2             0.49703    0.00864  57.55333   <1e-12
#>    v1~~v3             0.52388    0.00898  58.33147   <1e-12
#>    v1~~v4             0.48386    0.00841  57.53295   <1e-12
#>    v2                 1.00886    0.01427  70.71068         
#>    v3                 1.03278    0.01461  70.71068         
#>    v4                 0.95810    0.01355  70.71068         

## Using just sufficient statistics
n <- nrow(dd)
e0 <- estimate(m,data=list(S=cov(dd)*(n-1)/n,mu=colMeans(dd),n=n))
rm(dd)

## Multiple group analysis
m <- lvm()
regression(m) <- c(y1,y2,y3)~u
regression(m) <- u~x
d1 <- sim(m,100,p=c("u,u"=1,"u~x"=1))
d2 <- sim(m,100,p=c("u,u"=2,"u~x"=-1))

mm <- baptize(m)
regression(mm,u~x) <- NA
covariance(mm,~u) <- NA
intercept(mm,~u) <- NA
ee <- estimate(list(mm,mm),list(d1,d2))

## Missing data
d0 <- makemissing(d1,cols=1:2)
e0 <- estimate(m,d0,missing=TRUE)
e0
#>                     Estimate Std. Error  Z value Pr(>|z|)
#> Regressions:                                             
#>    y1~u              0.95061    0.07415 12.81979   <1e-12
#>     y2~u             1.11866    0.07109 15.73593   <1e-12
#>    y3~u              1.05029    0.06657 15.77822   <1e-12
#>     u~x              1.07240    0.10414 10.29807   <1e-12
#> Intercepts:                                              
#>    y1               -0.05948    0.11110 -0.53537   0.5924
#>    y2                0.03105    0.10751  0.28880   0.7727
#>    y3               -0.07114    0.10194 -0.69791   0.4852
#>    u                 0.02644    0.10703  0.24702   0.8049
#> Residual Variances:                                      
#>    y1                0.93568    0.14984  6.24435         
#>    y2                0.95506    0.14827  6.44139         
#>    y3                1.03629    0.14657  7.07041         
#>    u                 1.13503    0.16053  7.07042         
```
