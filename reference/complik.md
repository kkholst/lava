# Composite Likelihood for probit latent variable models

Estimate parameters in a probit latent variable model via a composite
likelihood decomposition.

## Usage

``` r
complik(
  x,
  data,
  k = 2,
  type = c("all", "nearest"),
  pairlist,
  messages = 0,
  estimator = "normal",
  quick = FALSE,
  ...
)
```

## Arguments

- x:

  `lvm`-object

- data:

  data.frame

- k:

  Size of composite groups

- type:

  Determines number of groups. With `type="nearest"` (default) only
  neighboring items will be grouped, e.g. for `k=2` (y1,y2),(y2,y3),...
  With `type="all"` all combinations of size `k` are included

- pairlist:

  A list of indices specifying the composite groups. Optional argument
  which overrides `k` and `type` but gives complete flexibility in the
  specification of the composite likelihood

- messages:

  Control amount of messages printed

- estimator:

  Model (pseudo-likelihood) to use for the pairs/groups

- quick:

  If TRUE the parameter estimates are calculated but all additional
  information such as standard errors are skipped

- ...:

  Additional arguments parsed on to lower-level functions

## Value

An object of class `estimate.complik` inheriting methods from `lvm`

## See also

estimate

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm(c(y1,y2,y3)~b*x+1*u[0],latent=~u)
ordinal(m,K=2) <- ~y1+y2+y3
d <- sim(m,50,seed=1)
if (requireNamespace("mets", quietly=TRUE)) {
   op <- lava.options(Dmethod="simple")
   e1 <- complik(m,d,control=list(trace=1),type="all")
   lava.options(op)
   e1
}
#>   0:     194.34186: 0.604480 0.584480 0.484480  0.00000 0.900000
#>   1:     163.15053: 0.540600 0.484448 0.209475 0.939309  1.06729
#>   2:     160.91018: 0.418339 0.334648 -0.188107  1.05711  1.14818
#>   3:     160.77829: 0.376957 0.300830 -0.155237  1.00706  1.18949
#>   4:     160.76097: 0.369541 0.291889 -0.159621  1.03751  1.20710
#>   5:     160.74923: 0.370653 0.293179 -0.167035  1.02015  1.23923
#>   6:     160.73021: 0.389084 0.308032 -0.176547  1.05300  1.30117
#>   7:     160.72405: 0.378467 0.292853 -0.153296  1.07323  1.36654
#>   8:     160.71872: 0.379911 0.296385 -0.160055  1.06751  1.36767
#>   9:     160.71668: 0.381244 0.299833 -0.167225  1.06316  1.36877
#>  10:     160.71655: 0.381339 0.300469 -0.169588  1.06363  1.36932
#>  11:     160.71647: 0.380400 0.300416 -0.170938  1.06481  1.37402
#>  12:     160.71645: 0.381360 0.299981 -0.170525  1.06481  1.37453
#>  13:     160.71644: 0.381181 0.300608 -0.170558  1.06511  1.37554
#>  14:     160.71644: 0.381382 0.300270 -0.170881  1.06528  1.37649
#>  15:     160.71644: 0.381626 0.300372 -0.170224  1.06543  1.37730
#>  16:     160.71643: 0.381497 0.300612 -0.170485  1.06571  1.37828
#>  17:     160.71643: 0.381626 0.300548 -0.170652  1.06574  1.37858
#>  18:     160.71643: 0.381584 0.300545 -0.170621  1.06576  1.37856
#>  19:     160.71643: 0.381588 0.300551 -0.170618  1.06575  1.37856
#>                     Estimate Std. Error  Z value  Pr(>|z|)
#> Regressions:                                              
#>    y1~x              1.06575    0.29344  3.63189 0.0002814
#> Intercepts:                                               
#>    y1                0.38159    0.31438  1.21377    0.2248
#>    y2                0.30055    0.29304  1.02563    0.3051
#>    y3               -0.17062    0.28519 -0.59825    0.5497
#> Residual Variances:                                       
#>    u                 1.37856    0.79666  1.73043          
```
