# Estimate probabilities in contingency table

Estimate probabilities in contingency table

## Usage

``` r
multinomial(
  x,
  data = parent.frame(),
  marginal = FALSE,
  transform,
  vcov = TRUE,
  IC = TRUE,
  ...
)
```

## Arguments

- x:

  Formula (or matrix or data.frame with observations, 1 or 2 columns)

- data:

  Optional data.frame

- marginal:

  If TRUE the marginals are estimated

- transform:

  Optional transformation of parameters (e.g., logit)

- vcov:

  Calculate asymptotic variance (default TRUE)

- IC:

  Return ic decomposition (default TRUE)

- ...:

  Additional arguments to lower-level functions

## Author

Klaus K. Holst

## Examples

``` r
set.seed(1)
breaks <- c(-Inf,-1,0,Inf)
m <- lvm(); covariance(m,pairwise=TRUE) <- ~y1+y2+y3+y4
d <- transform(sim(m,5e2),
              z1=cut(y1,breaks=breaks),
              z2=cut(y2,breaks=breaks),
              z3=cut(y3,breaks=breaks),
              z4=cut(y4,breaks=breaks))

multinomial(d[,5])
#> Call: multinomial(x = d[, 5])
#> 
#> Joint probabilities:
#> x
#> (-Inf,-1]    (-1,0]  (0, Inf] 
#>     0.154     0.350     0.496 
#> 
#>    Estimate Std.Err   2.5%  97.5%    P-value
#> p1    0.154 0.01614 0.1224 0.1856  1.425e-21
#> p2    0.350 0.02133 0.3082 0.3918  1.669e-60
#> p3    0.496 0.02236 0.4522 0.5398 5.068e-109
(a1 <- multinomial(d[,5:6]))
#> Call: multinomial(x = d[, 5:6])
#> 
#> Joint probabilities:
#>            z2
#> z1          (-Inf,-1] (-1,0] (0, Inf]
#>   (-Inf,-1]     0.064  0.062    0.028
#>   (-1,0]        0.066  0.146    0.138
#>   (0, Inf]      0.040  0.154    0.302
#> 
#> Conditional probabilities:
#>            z2
#> z1           (-Inf,-1]     (-1,0]   (0, Inf]
#>   (-Inf,-1] 0.41558442 0.40259740 0.18181818
#>   (-1,0]    0.18857143 0.41714286 0.39428571
#>   (0, Inf]  0.08064516 0.31048387 0.60887097
#> 
#>     Estimate  Std.Err    2.5%   97.5%   P-value
#> p11    0.064 0.010946 0.04255 0.08545 5.004e-09
#> p21    0.066 0.011104 0.04424 0.08776 2.780e-09
#> p31    0.040 0.008764 0.02282 0.05718 5.010e-06
#> p12    0.062 0.010785 0.04086 0.08314 8.986e-09
#> p22    0.146 0.015791 0.11505 0.17695 2.340e-20
#> p32    0.154 0.016142 0.12236 0.18564 1.425e-21
#> p13    0.028 0.007378 0.01354 0.04246 1.475e-04
#> p23    0.138 0.015424 0.10777 0.16823 3.657e-19
#> p33    0.302 0.020533 0.26176 0.34224 5.707e-49
(K1 <- kappa(a1)) ## Cohen's kappa
#>       Estimate Std.Err  2.5% 97.5%   P-value
#> kappa   0.2065 0.03547 0.137 0.276 5.805e-09

K2 <- kappa(d[,7:8])
## Testing difference K1-K2:
estimate(merge(K1,K2,id=TRUE),diff)
#>         Estimate Std.Err     2.5%  97.5% P-value
#> kappa.1  0.05756 0.04779 -0.03611 0.1512  0.2284

estimate(merge(K1,K2,id=FALSE),diff) ## Wrong std.err ignoring dependence
#>         Estimate Std.Err     2.5%  97.5% P-value
#> kappa.1  0.05756 0.04997 -0.04037 0.1555  0.2493
sqrt(vcov(K1)+vcov(K2))
#>            kappa
#> kappa 0.04996804

## Average of the two kappas:
estimate(merge(K1,K2,id=TRUE),function(x) mean(x))
#>    Estimate Std.Err   2.5%  97.5%  P-value
#> p1   0.2353 0.02603 0.1843 0.2863 1.57e-19
estimate(merge(K1,K2,id=FALSE),function(x) mean(x)) ## Independence
#>    Estimate Std.Err   2.5%  97.5%  P-value
#> p1   0.2353 0.02498 0.1863 0.2842 4.64e-21
##'
## Goodman-Kruskal's gamma
m2 <- lvm(); covariance(m2) <- y1~y2
breaks1 <- c(-Inf,-1,0,Inf)
breaks2 <- c(-Inf,0,Inf)
d2 <- transform(sim(m2,5e2),
              z1=cut(y1,breaks=breaks1),
              z2=cut(y2,breaks=breaks2))

(g1 <- gkgamma(d2[,3:4]))
#> Call: gkgamma(x = d2[, 3:4])
#> ────────────────────────────────────────────────────────────────────────────────
#> n = 500
#> 
#>       Estimate  Std.Err    2.5%   97.5%   P-value
#> C      0.26654 0.013898 0.23931 0.29378 5.522e-82
#> D      0.06619 0.007974 0.05056 0.08182 1.033e-16
#> gamma  0.60214 0.053796 0.49670 0.70757 4.411e-29
## same as
if (FALSE) { # \dontrun{
gkgamma(table(d2[,3:4]))
gkgamma(multinomial(d2[,3:4]))
} # }

##partial gamma
d2$x <- rbinom(nrow(d2),2,0.5)
gkgamma(z1~z2|x,data=d2)
#> Call: gkgamma(x = z1 ~ z2 | x, data = d2)
#> ────────────────────────────────────────────────────────────────────────────────
#> Strata:
#> 
#> 0 (n=112):
#>   Estimate Std.Err    2.5%   97.5%   P-value
#> C  0.29464 0.02999 0.23587 0.35342 8.758e-23
#> D  0.04624 0.01379 0.01921 0.07326 7.981e-04
#> 
#> 1 (n=248):
#>   Estimate Std.Err    2.5%  97.5%   P-value
#> C  0.24340 0.01958 0.20502 0.2818 1.797e-35
#> D  0.08126 0.01256 0.05665 0.1059 9.786e-11
#> 
#> 2 (n=140):
#>   Estimate Std.Err    2.5%  97.5%   P-value
#> C  0.28520 0.02604 0.23417 0.3362 6.415e-28
#> D  0.05806 0.01436 0.02992 0.0862 5.261e-05
#> 
#> ────────────────────────────────────────────────────────────────────────────────
#> 
#> n = 500
#> 
#> Gamma coefficient:
#> 
#>        Estimate Std.Err   2.5%  97.5%   P-value
#> γ:0      0.7287 0.09062 0.5511 0.9063 8.894e-16
#> γ:1      0.4994 0.08649 0.3299 0.6689 7.744e-09
#> γ:2      0.6617 0.09293 0.4796 0.8439 1.076e-12
#> pgamma   0.5663 0.06074 0.4473 0.6854 1.126e-20
```
