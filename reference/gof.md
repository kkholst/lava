# Extract model summaries and GOF statistics for model object

Calculates various GOF statistics for model object including global
chi-squared test statistic and AIC. Extract model-specific mean and
variance structure, residuals and various predicitions.

## Usage

``` r
gof(object, ...)

# S3 method for class 'lvmfit'
gof(object, chisq=FALSE, level=0.90, rmsea.threshold=0.05,all=FALSE,...)

moments(x,...)

# S3 method for class 'lvm'
moments(x, p, debug=FALSE, conditional=FALSE, data=NULL, latent=FALSE, ...)

# S3 method for class 'lvmfit'
logLik(object, p=coef(object),
                      data=model.frame(object),
                      model=object$estimator,
                      weights=Weights(object),
                      data2=object$data$data2,
                          ...)

# S3 method for class 'lvmfit'
score(x, data=model.frame(x), p=pars(x), model=x$estimator,
                   weights=Weights(x), data2=x$data$data2, ...)

# S3 method for class 'lvmfit'
information(x,p=pars(x),n=x$data$n,data=model.frame(x),
                   model=x$estimator,weights=Weights(x), data2=x$data$data2, ...)
```

## Arguments

- object:

  Model object

- ...:

  Additional arguments to be passed to the low level functions

- x:

  Model object

- p:

  Parameter vector used to calculate statistics

- data:

  Data.frame to use

- latent:

  If TRUE predictions of latent variables are included in output

- data2:

  Optional second data.frame (only for censored observations)

- weights:

  Optional weight matrix

- n:

  Number of observations

- conditional:

  If TRUE the conditional moments given the covariates are calculated.
  Otherwise the joint moments are calculated

- model:

  String defining estimator, e.g. "gaussian" (see `estimate`)

- debug:

  Debugging only

- chisq:

  Boolean indicating whether to calculate chi-squared goodness-of-fit
  (always TRUE for estimator='gaussian')

- level:

  Level of confidence limits for RMSEA

- rmsea.threshold:

  Which probability to calculate, Pr(RMSEA\<rmsea.treshold)

- all:

  Calculate all (ad hoc) FIT indices: TLI, CFI, NFI, SRMR, ...

## Value

A `htest`-object.

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm(list(y~v1+v2+v3+v4,c(v1,v2,v3,v4)~x))
set.seed(1)
dd <- sim(m,1000)
e <- estimate(m, dd)
gof(e,all=TRUE,rmsea.threshold=0.05,level=0.9)
#> 
#>  Number of observations = 1000 
#>  BIC = 14585.57 
#>  AIC = 14468.26 
#>  log-Likelihood of model = -7216.128 
#> 
#>  log-Likelihood of saturated model = -7212.5 
#>  Chi-squared statistic: q = 7.254653 , df = 7 
#>   P(Q>q) = 0.4028559 
#> 
#>  RMSEA (90% CI): 0.006 (0;0.0397)
#>   P(RMSEA<0.05)=0.9916145
#>  TLI = 0.9998998 
#>  CFI = 0.9999532 
#>  NFI = 0.9986715 
#>  SRMR = 0.008682085 
#> 
#> rank(Information) = 18 (p=18)
#> condition(Information) = 10.37525
#> mean(score^2) = 4.216767e-09 


set.seed(1)
m <- lvm(list(c(y1,y2,y3)~u,y1~x)); latent(m) <- ~u
regression(m,c(y2,y3)~u) <- "b"
d <- sim(m,1000)
e <- estimate(m,d)
rsq(e)
#> $`R-squared`
#>           y1           y2           y3            u 
#> 6.714238e-01 5.109812e-01 5.276472e-01 1.443290e-15 
#> 
#> $`Variance explained by 'u'`
#>        y1        y2        y3 
#> 0.3697894 0.5109812 0.5276472 
#> 
##'
rr <- rsq(e,TRUE)
rr
#> 
#> R-squared:
#> 
#>     Estimate    Std.Err      2.5%     97.5%       P-value
#> y1 0.6666507 0.02449714 0.6186372 0.7146642 4.506818e-163
#> y2 0.5062724 0.02751655 0.4523409 0.5602038  1.342309e-75
#> y3 0.5319590 0.02627482 0.4804613 0.5834567  3.855758e-91
estimate(rr,contrast=rbind(c(1,-1,0),c(1,0,-1),c(0,1,-1)))
#>             Estimate Std.Err     2.5%   97.5%   P-value
#> [y1] - [y2]  0.16038 0.04040  0.08119 0.23956 7.197e-05
#> [y1] - [y3]  0.13469 0.03884  0.05857 0.21081 5.244e-04
#> [y2] - [y3] -0.02569 0.02786 -0.08029 0.02891 3.565e-01
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [y1] - [y2] = 0
#>   [y1] - [y3] = 0
#>   [y2] - [y3] = 0 
#>  
#> chisq = 16.2844, df = 2, p-value = 0.000291
```
