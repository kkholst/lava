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
#>    y~v1               1.00125    0.01784  56.10838   <1e-12
#>    y~v2               0.99412    0.01085  91.63608   <1e-12
#>    y~v3               1.00628    0.01103  91.22031   <1e-12
#>    y~v4               0.99382    0.01099  90.44727   <1e-12
#>     v1~x              0.99191    0.01008  98.39985   <1e-12
#>    v2~x               1.00066    0.01011  98.97991   <1e-12
#>     v3~x              0.98852    0.01020  96.95461   <1e-12
#>    v4~x               1.00548    0.00986 101.97506   <1e-12
#> Intercepts:                                                
#>    y                  0.00592    0.01001   0.59091   0.5546
#>    v1                -0.00891    0.01001  -0.88995   0.3735
#>    v2                 0.00576    0.01004   0.57354   0.5663
#>    v3                -0.01314    0.01013  -1.29775   0.1944
#>    v4                -0.01015    0.00979  -1.03664   0.2999
#> Residual Variances:                                        
#>    y                  1.00179    0.01417  70.71068         
#>    v1                 1.00266    0.01121  89.48147         
#>    v1~~v2             0.49740    0.00864  57.55765   <1e-12
#>    v1~~v3             0.52049    0.00893  58.26125   <1e-12
#>    v1~~v4             0.48318    0.00841  57.48012   <1e-12
#>    v2                 1.00850    0.01426  70.71068         
#>    v3                 1.02572    0.01451  70.71068         
#>    v4                 0.95929    0.01357  70.71068         

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
#>    y1~u              1.07277    0.07428 14.44300   <1e-12
#>     y2~u             1.03533    0.07713 13.42284   <1e-12
#>    y3~u              1.04238    0.06356 16.39954   <1e-12
#>     u~x              1.12259    0.12611  8.90200   <1e-12
#> Intercepts:                                              
#>    y1               -0.02455    0.11003 -0.22312   0.8234
#>    y2               -0.02719    0.10900 -0.24948    0.803
#>    y3                0.06997    0.09110  0.76812   0.4424
#>    u                -0.04358    0.10698 -0.40737   0.6837
#> Residual Variances:                                      
#>    y1                1.01596    0.15678  6.48011         
#>    y2                0.95015    0.15025  6.32389         
#>    y3                0.82857    0.11719  7.07021         
#>    u                 1.14417    0.16183  7.07044         
```
