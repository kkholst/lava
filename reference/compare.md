# Statistical tests

Performs Likelihood-ratio, Wald and score tests

## Usage

``` r
compare(object, ...)
```

## Arguments

- object:

  `lvmfit`-object

- ...:

  Additional arguments to low-level functions

## Value

Matrix of test-statistics and p-values

## See also

[`modelsearch`](http://kkholst.github.io/lava/reference/modelsearch.md),
[`equivalence`](http://kkholst.github.io/lava/reference/equivalence.md)

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm();
regression(m) <- c(y1,y2,y3) ~ eta; latent(m) <- ~eta
regression(m) <- eta ~ x
m2 <- regression(m, c(y3,eta) ~ x)
set.seed(1)
d <- sim(m,1000)
e <- estimate(m,d)
e2 <- estimate(m2,d)

compare(e)
#> 
#>  - Likelihood ratio test -
#> 
#> data:  
#> chisq = 2.7373, df = 2, p-value = 0.2544
#> sample estimates:
#>           log likelihood (model) log likelihood (saturated model) 
#>                        -5045.863                        -5044.494 
#> 

compare(e,e2) ## LRT, H0: y3<-x=0
#> 
#>  - Likelihood ratio test -
#> 
#> data:  
#> chisq = 1.6297, df = 1, p-value = 0.2017
#> sample estimates:
#> log likelihood (model 1) log likelihood (model 2) 
#>                -5045.863                -5045.048 
#> 
compare(e,scoretest=y3~x) ## Score-test, H0: y3~x=0
#> 
#>  - Score test -
#> 
#> data:  y3 ~ x
#> chisq = 1.6059, df = 1, p-value = 0.2051
#> 
compare(e2,par=c("y3~x")) ## Wald-test, H0: y3~x=0
#> 
#>  - Wald test -
#> 
#>  Null Hypothesis:
#>  [y3~x] = 0
#> 
#> data:  
#> chisq = 1.5752, df = 1, p-value = 0.2095
#> sample estimates:
#>           Estimate    Std.Err     2.5%      97.5%
#> [y3~x] -0.08157255 0.06499477 -0.20896 0.04581487
#> 

B <- diag(2); colnames(B) <- c("y2~eta","y3~eta")
compare(e2,contrast=B,null=c(1,1))
#> 
#>  - Wald test -
#> 
#>  Null Hypothesis:
#>  [y2~eta] = 1
#>  [y3~eta] = 1
#> 
#> data:  
#> chisq = 0.40264, df = 2, p-value = 0.8177
#> sample estimates:
#>          Estimate    Std.Err      2.5%    97.5%
#> [y2~eta] 1.019845 0.03770718 0.9459398 1.093749
#> [y3~eta] 1.028685 0.05598807 0.9189509 1.138420
#> 

B <- rep(0,length(coef(e2))); B[1:3] <- 1
compare(e2,contrast=B)
#> 
#>  - Wald test -
#> 
#>  Null Hypothesis:
#>  [y2] + [y3] + [eta] = 0
#> 
#> data:  
#> chisq = 0.15653, df = 1, p-value = 0.6924
#> sample estimates:
#>                       Estimate    Std.Err       2.5%     97.5%
#> [y2] + [y3] + [eta] 0.02605068 0.06584406 -0.1030013 0.1551027
#> 

compare(e,scoretest=list(y3~x,y2~x))
#> 
#>  - Score test -
#> 
#> data:  y3 ~ xy2 ~ x
#> chisq = 2.7607, df = 2, p-value = 0.2515
#> 
```
