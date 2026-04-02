# Calculate confidence limits for parameters

Calculate Wald og Likelihood based (profile likelihood) confidence
intervals

## Usage

``` r
# S3 method for class 'lvmfit'
confint(
  object,
  parm = seq_len(length(coef(object))),
  level = 0.95,
  profile = FALSE,
  curve = FALSE,
  n = 20,
  interval = NULL,
  lower = TRUE,
  upper = TRUE,
  ...
)
```

## Arguments

- object:

  `lvm`-object.

- parm:

  Index of which parameters to calculate confidence limits for.

- level:

  Confidence level

- profile:

  Logical expression defining whether to calculate confidence limits via
  the profile log likelihood

- curve:

  if FALSE and profile is TRUE, confidence limits are returned.
  Otherwise, the profile curve is returned.

- n:

  Number of points to evaluate profile log-likelihood in over the
  interval defined by `interval`

- interval:

  Interval over which the profiling is done

- lower:

  If FALSE the lower limit will not be estimated (profile intervals
  only)

- upper:

  If FALSE the upper limit will not be estimated (profile intervals
  only)

- ...:

  Additional arguments to be passed to the low level functions

## Value

A 2xp matrix with columns of lower and upper confidence limits

## Details

Calculates either Wald confidence limits: \$\$\hat{\theta} \pm
z\_{\alpha/2}\*\hat\sigma\_{\hat\theta}\$\$ or profile likelihood
confidence limits, defined as the set of value \\\tau\\:
\$\$logLik(\hat\theta\_{\tau},\tau)-logLik(\hat\theta)\<
q\_{\alpha}/2\$\$

where \\q\_{\alpha}\\ is the \\\alpha\\ fractile of the \\\chi^2_1\\
distribution, and \\\hat\theta\_{\tau}\\ are obtained by maximizing the
log-likelihood with tau being fixed.

## See also

[`bootstrap`](https://kkholst.github.io/lava/reference/bootstrap.md)`{lvm}`

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm(y~x)
d <- sim(m,100)
e <- estimate(lvm(y~x), d)
confint(e,3,profile=TRUE)
#>          2.5 %   97.5 %
#> y~~y 0.5888514 1.026329
confint(e,3)
#>          2.5 %    97.5 %
#> y~~y 0.5547589 0.9802276
 ## Reduce Ex.timings
B <- bootstrap(e,R=50)
B
#> Non-parametric bootstrap statistics (R=50):
#> 
#>      Estimate     Bias         Std.Err      2.5 %        97.5 %      
#> y     0.124504326 -0.007753062  0.073936733  0.008375347  0.288089457
#> y~x   1.061750034  0.013940880  0.074844578  0.949601761  1.235976345
#> y~~y  0.767493210 -0.026454597  0.100031344  0.584892411  0.986727053
#> 
```
