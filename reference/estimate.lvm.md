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
  [`help(lava.options)`](http://kkholst.github.io/lava/reference/lava.options.md))

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
#>    y~v1               1.00086    0.01785  56.06986   <1e-12
#>    y~v2               0.99394    0.01085  91.59488   <1e-12
#>    y~v3               1.00627    0.01103  91.20677   <1e-12
#>    y~v4               0.99443    0.01099  90.47824   <1e-12
#>     v1~x              0.99205    0.01009  98.33296   <1e-12
#>    v2~x               1.00074    0.01012  98.90608   <1e-12
#>     v3~x              0.98866    0.01021  96.87999   <1e-12
#>    v4~x               1.00562    0.00987 101.88095   <1e-12
#> Intercepts:                                                
#>    y                  0.00609    0.01001   0.60833    0.543
#>    v1                -0.00912    0.01001  -0.91110   0.3622
#>    v2                 0.00566    0.01004   0.56365    0.573
#>    v3                -0.01290    0.01013  -1.27400   0.2027
#>    v4                -0.01082    0.00980  -1.10458   0.2693
#> Residual Variances:                                        
#>    y                  1.00204    0.01417  70.71068         
#>    v1                 1.00257    0.01120  89.48388         
#>    v1~~v2             0.49732    0.00864  57.55755   <1e-12
#>    v1~~v3             0.52040    0.00893  58.25942   <1e-12
#>    v1~~v4             0.48340    0.00841  57.48707   <1e-12
#>    v2                 1.00842    0.01426  70.71068         
#>    v3                 1.02582    0.01451  70.71068         
#>    v4                 0.95968    0.01357  70.71068         

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
#>    y1~u              1.08224    0.08055 13.43642   <1e-12
#>     y2~u             0.94896    0.06791 13.97290   <1e-12
#>    y3~u              1.07444    0.06228 17.25066   <1e-12
#>     u~x              1.12905    0.11653  9.68858   <1e-12
#> Intercepts:                                              
#>    y1               -0.08777    0.11215 -0.78267   0.4338
#>    y2                0.01866    0.10486  0.17799   0.8587
#>    y3                0.11980    0.09151  1.30921   0.1905
#>    u                -0.13072    0.10531 -1.24130   0.2145
#> Residual Variances:                                      
#>    y1                1.04880    0.16185  6.48015         
#>    y2                0.87854    0.13892  6.32385         
#>    y3                0.83338    0.11787  7.07023         
#>    u                 1.10811    0.15672  7.07044         
```
