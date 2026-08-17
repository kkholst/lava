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
   e1 <- complik(m,d,control=list(trace=1),type="all")
}
#>   0:     194.34186: 0.604480 0.584480 0.484480  0.00000 0.900000
#>   1:     163.15053: 0.540600 0.484448 0.209475 0.939309  1.06729
#>   2:     160.91018: 0.418339 0.334648 -0.188107  1.05711  1.14818
#>   3:     160.77829: 0.376957 0.300830 -0.155237  1.00706  1.18949
#>   4:     160.76097: 0.369541 0.291889 -0.159621  1.03751  1.20710
#>   5:     160.74923: 0.370653 0.293179 -0.167035  1.02015  1.23923
#>   6:     160.73532: 0.376191 0.296936 -0.168307  1.04656  1.26468
#>   7:     160.72558: 0.379780 0.299688 -0.168631  1.04345  1.30158
#>   8:     160.72091: 0.369539 0.303131 -0.165350  1.05346  1.33569
#>   9:     160.71897: 0.382701 0.295562 -0.170712  1.05926  1.33925
#>  10:     160.71733: 0.379230 0.296353 -0.167731  1.05972  1.35609
#>  11:     160.71701: 0.379746 0.299743 -0.170560  1.06152  1.35778
#>  12:     160.71683: 0.380300 0.298684 -0.169147  1.06113  1.36250
#>  13:     160.71663: 0.380886 0.299718 -0.170369  1.06381  1.36647
#>  14:     160.71650: 0.380928 0.300399 -0.170596  1.06425  1.37148
#>  15:     160.71645: 0.381806 0.299785 -0.170466  1.06543  1.37630
#>  16:     160.71643: 0.381561 0.300967 -0.170583  1.06589  1.37944
#>  17:     160.71643: 0.381542 0.300509 -0.170678  1.06575  1.37858
#>  18:     160.71643: 0.381620 0.300565 -0.170579  1.06575  1.37852
#>  19:     160.71643: 0.381587 0.300550 -0.170619  1.06575  1.37856
#>  20:     160.71643: 0.381587 0.300550 -0.170619  1.06575  1.37856
```
