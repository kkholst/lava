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
#>    y~v1               1.00124    0.01784  56.11097   <1e-12
#>    y~v2               0.99419    0.01085  91.63349   <1e-12
#>    y~v3               1.00648    0.01103  91.24373   <1e-12
#>    y~v4               0.99339    0.01099  90.41890   <1e-12
#>     v1~x              0.99124    0.01008  98.32354   <1e-12
#>    v2~x               1.00029    0.01011  98.92904   <1e-12
#>     v3~x              0.98819    0.01020  96.91956   <1e-12
#>    v4~x               1.00491    0.00986 101.89550   <1e-12
#> Intercepts:                                                
#>    y                  0.00589    0.01001   0.58805   0.5565
#>    v1                -0.00886    0.01001  -0.88508   0.3761
#>    v2                 0.00558    0.01004   0.55588   0.5783
#>    v3                -0.01287    0.01013  -1.27047   0.2039
#>    v4                -0.01026    0.00979  -1.04729    0.295
#> Residual Variances:                                        
#>    y                  1.00171    0.01417  70.71068         
#>    v1                 1.00251    0.01120  89.47976         
#>    v1~~v2             0.49749    0.00864  57.56218   <1e-12
#>    v1~~v3             0.52029    0.00893  58.25756   <1e-12
#>    v1~~v4             0.48308    0.00840  57.47552   <1e-12
#>    v2                 1.00843    0.01426  70.71068         
#>    v3                 1.02543    0.01450  70.71068         
#>    v4                 0.95938    0.01357  70.71068         

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
#>    y1~u              1.06614    0.07743 13.76980   <1e-12
#>     y2~u             1.00643    0.07091 14.19311   <1e-12
#>    y3~u              1.01669    0.06287 16.17022   <1e-12
#>     u~x              1.17785    0.11892  9.90456   <1e-12
#> Intercepts:                                              
#>    y1               -0.03623    0.11353 -0.31908   0.7497
#>    y2                0.02718    0.11168  0.24341   0.8077
#>    y3                0.06231    0.09336  0.66737   0.5045
#>    u                -0.06776    0.10529 -0.64356   0.5199
#> Residual Variances:                                      
#>    y1                1.07366    0.16569  6.48014         
#>    y2                0.97608    0.15532  6.28423         
#>    y3                0.86575    0.12245  7.07022         
#>    u                 1.10552    0.15636  7.07040         
```
