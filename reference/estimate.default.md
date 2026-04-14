# Estimation of functional of parameters

Estimation of functional of parameters. Wald tests, robust standard
errors, cluster robust standard errors

## Usage

``` r
# Default S3 method
estimate(
  x = NULL,
  f = NULL,
  ...,
  data,
  id,
  stack = TRUE,
  average = FALSE,
  subset,
  score.deriv,
  level = 0.95,
  IC = robust,
  type = c("robust", "df", "mbn"),
  var.adj,
  keep,
  use,
  regex = FALSE,
  ignore.case = FALSE,
  contrast,
  null,
  vcov,
  coef,
  robust = TRUE,
  df = NULL,
  print = NULL,
  labels,
  label.width,
  only.coef = FALSE,
  back.transform = NULL,
  folds = 0,
  R = 0,
  null.sim
)
```

## Arguments

- x:

  model object (`glm`, `lvmfit`, ...)

- f:

  transformation of model parameters and (optionally) data, or contrast
  matrix (or vector)

- ...:

  additional arguments to lower level functions

- data:

  `data.frame`

- id:

  (optional) id-variable corresponding to ic decomposition of model
  parameters.

- stack:

  if TRUE (default) the i.i.d. decomposition is automatically stacked
  according to 'id'

- average:

  if TRUE averages are calculated

- subset:

  (optional) subset of data.frame on which to condition (logical
  expression or variable name)

- score.deriv:

  (optional) derivative of mean score function

- level:

  level of confidence limits

- IC:

  if TRUE (default) the influence function decompositions are also
  returned (extract with `IC` method)

- type:

  type of small-sample correction

- var.adj:

  variance adjustment parameter for small-sample correction

- keep:

  (optional) index of parameters to keep from final result

- use:

  (optional) index of parameters to use in calculations

- regex:

  If TRUE use regular expression (perl compatible) for keep, use
  arguments

- ignore.case:

  Ignore case-sensitiveness in regular expression

- contrast:

  (optional) Contrast matrix for final Wald test

- null:

  (optional) null hypothesis to test

- vcov:

  (optional) covariance matrix of parameter estimates

- coef:

  (optional) parameter coefficient

- robust:

  if TRUE robust standard errors are calculated

- df:

  degrees of freedom (default obtained from 'df.residual')

- print:

  (optional) print function for the resulting estimate objecta

- labels:

  (optional) names of coefficients

- label.width:

  (optional) max width of labels

- only.coef:

  if TRUE only the coefficient matrix is return

- back.transform:

  (optional) transform of parameters and confidence intervals

- folds:

  (optional) aggregate influence functions (divide and conquer)

- R:

  Number of simulations (simulated p-values)

- null.sim:

  Mean under the null for simulations

## Details

influence function decomposition of estimator \\\widehat{\theta}\\ based
on data \\Z_1,\ldots,Z_n\\: \$\$\sqrt{n}(\widehat{\theta}-\theta) =
\frac{1}{\sqrt{n}}\sum\_{i=1}^n IC(Z_i; P) + o_p(1)\$\$ can be extracted
with the `IC` method.

## See also

estimate.array

## Examples

``` r
## Simulation from logistic regression model
m <- lvm(y~x+z);
distribution(m,y~x) <- binomial.lvm("logit")
d <- sim(m,1000)
g <- glm(y~z+x,data=d,family=binomial())
g0 <- glm(y~1,data=d,family=binomial())

## LRT
estimate(g, g0)
#> 
#>  - Likelihood ratio test -
#> 
#> data:  
#> chisq = 209.91, df = 2, p-value < 2.2e-16
#> sample estimates:
#> log likelihood (model 1) log likelihood (model 2) 
#>                -567.6493                -672.6041 
#> 

## Plain estimates (robust standard errors)
estimate(g)
#>              Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept) -0.001888 0.09785 -0.1937 0.1899 9.846e-01
#> z            0.953974 0.08114  0.7949 1.1130 6.481e-32
#> x            1.009058 0.14720  0.7205 1.2976 7.135e-12

## Testing contrasts
estimate(g, null=0)
#>              Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept) -0.001888 0.09785 -0.1937 0.1899 9.846e-01
#> z            0.953974 0.08114  0.7949 1.1130 6.481e-32
#> x            1.009058 0.14720  0.7205 1.2976 7.135e-12
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [(Intercept)] = 0
#>   [z] = 0
#>   [x] = 0 
#>  
#> chisq = 180.732, df = 3, p-value < 2.2e-16
estimate(g, rbind(c(1,1,0), c(1,0,2)))
#>                      Estimate Std.Err   2.5% 97.5%   P-value
#> [(Intercept)] + [z]    0.9521  0.1260 0.7052 1.199 4.075e-14
#> [(Intercept)] + 2[x]   2.0162  0.2404 1.5451 2.487 4.939e-17
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [(Intercept)] + [z] = 0
#>   [(Intercept)] + 2[x] = 0 
#>  
#> chisq = 153.6241, df = 2, p-value < 2.2e-16
estimate(g, rbind(c(1,1,0), c(1,0,2)), null=c(1,2))
#>                      Estimate Std.Err   2.5% 97.5% P-value
#> [(Intercept)] + [z]    0.9521  0.1260 0.7052 1.199  0.7037
#> [(Intercept)] + 2[x]   2.0162  0.2404 1.5451 2.487  0.9462
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [(Intercept)] + [z] = 1
#>   [(Intercept)] + 2[x] = 2 
#>  
#> chisq = 0.1447, df = 2, p-value = 0.9302
estimate(g, 2:3) ## same as cbind(0,1,-1)
#>           Estimate Std.Err    2.5%  97.5% P-value
#> [z] - [x] -0.05508  0.1537 -0.3563 0.2461    0.72
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [z] - [x] = 0 
#>  
#> chisq = 0.1285, df = 1, p-value = 0.72
estimate(g, as.list(2:3)) ## same as rbind(c(0,1,0),c(0,0,1))
#>   Estimate Std.Err   2.5% 97.5%   P-value
#> z    0.954 0.08114 0.7949 1.113 6.481e-32
#> x    1.009 0.14720 0.7205 1.298 7.135e-12
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [z] = 0
#>   [x] = 0 
#>  
#> chisq = 159.9584, df = 2, p-value < 2.2e-16
## Alternative syntax
estimate(g, "z", "z"-"x", 2*"z"-3*"x")
#>             Estimate Std.Err    2.5%   97.5%   P-value
#> z            0.95397 0.08114  0.7949  1.1130 6.481e-32
#> [z] - [x]   -0.05508 0.15367 -0.3563  0.2461 7.200e-01
#> 2[z] - 3[x] -1.11922 0.43991 -1.9814 -0.2570 1.095e-02
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [z] = 0
#>   [z] - [x] = 0
#>   2[z] - 3[x] = 0 
#>  
#> chisq = 159.9584, df = 2, p-value < 2.2e-16
estimate(g, "?")  ## Wildcards
#>           Estimate Std.Err    2.5%  97.5% P-value
#> [z] - [x] -0.05508  0.1537 -0.3563 0.2461    0.72
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [z] - [x] = 0 
#>  
#> chisq = 0.1285, df = 1, p-value = 0.72
estimate(g, "*Int*", "z")
#>              Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept) -0.001888 0.09785 -0.1937 0.1899 9.846e-01
#> z            0.953974 0.08114  0.7949 1.1130 6.481e-32
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [(Intercept)] = 0
#>   [z] = 0 
#>  
#> chisq = 138.2721, df = 2, p-value < 2.2e-16
estimate(g, "1", "2"-"3", null = c(0,1))
#>              Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept) -0.001888 0.09785 -0.1937 0.1899 9.846e-01
#> [z] - [x]   -0.055083 0.15367 -0.3563 0.2461 6.606e-12
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [(Intercept)] = 0
#>   [z] - [x] = 1 
#>  
#> chisq = 77.8686, df = 2, p-value < 2.2e-16
estimate(g, 2, 3)
#>   Estimate Std.Err   2.5% 97.5%   P-value
#> z    0.954 0.08114 0.7949 1.113 6.481e-32
#> x    1.009 0.14720 0.7205 1.298 7.135e-12
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [z] = 0
#>   [x] = 0 
#>  
#> chisq = 159.9584, df = 2, p-value < 2.2e-16

## Usual (non-robust) confidence intervals
estimate(g, robust=FALSE)
#>              Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept) -0.001888 0.09817 -0.1943 0.1905 9.847e-01
#> z            0.953974 0.08318  0.7909 1.1170 1.892e-30
#> x            1.009058 0.14639  0.7221 1.2960 5.469e-12

## Transformations
estimate(g, function(p) p[1]+p[2])
#>             Estimate Std.Err   2.5% 97.5%   P-value
#> (Intercept)   0.9521   0.126 0.7052 1.199 4.075e-14

## Multiple parameters
e <- estimate(g, function(p) c(p[1]+p[2], p[1]*p[2]))
e
#>                Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept)    0.952086 0.12596  0.7052 1.1990 4.075e-14
#> (Intercept).1 -0.001801 0.09335 -0.1848 0.1812 9.846e-01
vcov(e)
#>             (Intercept) (Intercept)
#> (Intercept)  0.01586612 0.008982930
#> (Intercept)  0.00898293 0.008714966

## Label new parameters
estimate(g, function(p) list("a1"=p[1]+p[2], "b1"=p[1]*p[2]))
#>     Estimate Std.Err    2.5%  97.5%   P-value
#> a1  0.952086 0.12596  0.7052 1.1990 4.075e-14
#> b1 -0.001801 0.09335 -0.1848 0.1812 9.846e-01
##'
## Multiple group
m <- lvm(y~x)
m <- baptize(m)
d2 <- d1 <- sim(m,50,seed=1)
e <- estimate(list(m,m),list(d1,d2))
estimate(e) ## Wrong
#>        Estimate Std.Err     2.5%  97.5%   P-value
#> y@1      0.1044 0.08277 -0.05785 0.2666 2.073e-01
#> y~x@1    0.9665 0.08727  0.79541 1.1375 1.677e-28
#> y~~y@1   0.6764 0.10629  0.46803 0.8847 1.977e-10
ee <- estimate(e, id=rep(seq(nrow(d1)), 2)) ## Clustered
ee
#>        Estimate Std.Err    2.5%  97.5%   P-value
#> y@1      0.1044  0.1171 -0.1250 0.3338 3.725e-01
#> y~x@1    0.9665  0.1234  0.7246 1.2084 4.859e-15
#> y~~y@1   0.6764  0.1503  0.3817 0.9710 6.814e-06
estimate(lm(y~x,d1))
#>             Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept)   0.1044  0.1171 -0.1251 0.3338 3.726e-01
#> x             0.9665  0.1234  0.7246 1.2084 4.853e-15

## Marginalize
f <- function(p,data)
  list(p0=lava:::expit(p["(Intercept)"] + p["z"]*data[,"z"]),
       p1=lava:::expit(p["(Intercept)"] + p["x"] + p["z"]*data[,"z"]))
e <- estimate(g, f, average=TRUE)
e
#>    Estimate Std.Err   2.5%  97.5%    P-value
#> p0   0.5010 0.02140 0.4591 0.5429 3.025e-121
#> p1   0.7007 0.01973 0.6620 0.7393 3.516e-276
estimate(e,diff)
#>    Estimate Std.Err  2.5%  97.5%   P-value
#> p1   0.1997  0.0279 0.145 0.2543 8.194e-13
estimate(e,cbind(1,1))
#>             Estimate Std.Err  2.5% 97.5% P-value
#> [p0] + [p1]    1.202 0.03027 1.142 1.261       0
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [p0] + [p1] = 0 
#>  
#> chisq = 1576.103, df = 1, p-value < 2.2e-16

## Clusters and subset (conditional marginal effects)
d$id <- rep(seq(nrow(d)/4),each=4)
estimate(g,function(p,data)
         list(p0=lava:::expit(p[1] + p["z"]*data[,"z"])),
         subset=d$z>0, id=d$id, average=TRUE)
#>    Estimate Std.Err   2.5%  97.5%    P-value
#> p0   0.6754 0.02282 0.6307 0.7202 1.558e-192

## More examples with clusters:
m <- lvm(c(y1,y2,y3)~u+x)
d <- sim(m,10)
l1 <- glm(y1~x,data=d)
l2 <- glm(y2~x,data=d)
l3 <- glm(y3~x,data=d)

## Some random id-numbers
id1 <- c(1,1,4,1,3,1,2,3,4,5)
id2 <- c(1,2,3,4,5,6,7,8,1,1)
id3 <- seq(10)

## Un-stacked and stacked i.i.d. decomposition
IC(estimate(l1,id=id1,stack=FALSE))
#>   (Intercept)           x
#> 1 -0.04904705 -0.14189362
#> 1  0.14744213 -0.15724052
#> 4  1.15168804 -0.54934034
#> 1  0.37775868  0.06124426
#> 3 -0.01647188 -0.03981494
#> 1 -1.75658710  2.38806661
#> 2 -0.48187968 -0.18690182
#> 3 -0.56660619 -0.30858145
#> 4  0.77907752 -0.82047631
#> 5  0.41462553 -0.24506186
#> attr(,"bread")
#>             (Intercept)          x
#> (Intercept)   0.6136517 -0.1000796
#> x            -0.1000796  0.9011828
IC(estimate(l1,id=id1))
#>   (Intercept)           x
#> 1  -0.6402167  1.07508836
#> 2  -0.2409398 -0.09345091
#> 3  -0.2915390 -0.17419820
#> 4   0.9653828 -0.68490833
#> 5   0.2073128 -0.12253093
#> attr(,"bread")
#>             (Intercept)          x
#> (Intercept)   0.6136517 -0.1000796
#> x            -0.1000796  0.9011828
#> attr(,"N")
#> [1] 10

## Combined i.i.d. decomposition
e1 <- estimate(l1,id=id1)
e2 <- estimate(l2,id=id2)
e3 <- estimate(l3,id=id3)
(a2 <- merge(e1,e2,e3))
#>               Estimate Std.Err     2.5% 97.5%   P-value
#> (Intercept)     0.5418  0.2472  0.05731 1.026 0.0283949
#> x               0.9636  0.2592  0.45568 1.472 0.0002006
#> (Intercept).1   0.2156  0.5105 -0.78505 1.216 0.6728124
#> x.1             1.2369  0.4854  0.28554 2.188 0.0108277
#> (Intercept).2   0.4769  0.2979 -0.10705 1.061 0.1094554
#> x.2             1.1545  0.3714  0.42660 1.882 0.0018795

## If all models were estimated on the same data we could use the
## syntax:
## Reduce(merge,estimate(list(l1,l2,l3)))

## Same:
IC(a1 <- merge(l1,l2,l3,id=list(id1,id2,id3)))
#>    (Intercept)          x (Intercept).1         x.1 (Intercept).2         x.2
#> 1   -1.2804333  2.1501767    3.53262530 -3.15431642     0.7554987  2.18566521
#> 2   -0.4818797 -0.1869018    0.06810688 -0.07263298     0.5200220 -0.55458052
#> 3   -0.5830781 -0.3483964    1.08830429 -0.51910711    -0.1215704  0.05798750
#> 4    1.9307656 -1.3698167   -2.27428000 -0.36871845    -1.5561326 -0.25228854
#> 5    0.4146255 -0.2450619    0.25612485  0.61909106    -0.4880168 -1.17960761
#> 6    0.0000000  0.0000000   -2.62738927  3.57191544     1.6742259 -2.27609724
#> 7    0.0000000  0.0000000    0.33520853  0.13001396     0.8010558  0.31069745
#> 8    0.0000000  0.0000000   -0.37870058 -0.20624550     0.1460747  0.07955428
#> 9    0.0000000  0.0000000    0.00000000  0.00000000    -1.3102941  1.37992084
#> 10   0.0000000  0.0000000    0.00000000  0.00000000    -0.4208633  0.24874863

IC(merge(l1,l2,l3,id=TRUE)) # one-to-one (same clusters)
#>    (Intercept)           x (Intercept).1         x.1 (Intercept).2         x.2
#> 1  -0.04904705 -0.14189362   0.004626991  0.01338593     0.7554987  2.18566521
#> 2   0.14744213 -0.15724052   0.068106877 -0.07263298     0.5200220 -0.55458052
#> 3   1.15168804 -0.54934034   1.088304293 -0.51910711    -0.1215704  0.05798750
#> 4   0.37775868  0.06124426  -2.274280005 -0.36871845    -1.5561326 -0.25228854
#> 5  -0.01647188 -0.03981494   0.256124851  0.61909106    -0.4880168 -1.17960761
#> 6  -1.75658710  2.38806661  -2.627389265  3.57191544     1.6742259 -2.27609724
#> 7  -0.48187968 -0.18690182   0.335208530  0.13001396     0.8010558  0.31069745
#> 8  -0.56660619 -0.30858145  -0.378700581 -0.20624550     0.1460747  0.07955428
#> 9   0.77907752 -0.82047631   2.342596539 -2.46707796    -1.3102941  1.37992084
#> 10  0.41462553 -0.24506186   1.185401770 -0.70062440    -0.4208633  0.24874863
IC(merge(l1,l2,l3,id=FALSE)) # independence
#>    (Intercept)          x (Intercept).1         x.1 (Intercept).2        x.2
#> 1  -0.14714116 -0.4256809    0.00000000  0.00000000     0.0000000  0.0000000
#> 2   0.44232639 -0.4717216    0.00000000  0.00000000     0.0000000  0.0000000
#> 3   3.45506412 -1.6480210    0.00000000  0.00000000     0.0000000  0.0000000
#> 4   1.13327605  0.1837328    0.00000000  0.00000000     0.0000000  0.0000000
#> 5  -0.04941565 -0.1194448    0.00000000  0.00000000     0.0000000  0.0000000
#> 6  -5.26976131  7.1641998    0.00000000  0.00000000     0.0000000  0.0000000
#> 7  -1.44563904 -0.5607055    0.00000000  0.00000000     0.0000000  0.0000000
#> 8  -1.69981856 -0.9257444    0.00000000  0.00000000     0.0000000  0.0000000
#> 9   2.33723256 -2.4614289    0.00000000  0.00000000     0.0000000  0.0000000
#> 10  1.24387660 -0.7351856    0.00000000  0.00000000     0.0000000  0.0000000
#> 11  0.00000000  0.0000000    0.01388097  0.04015779     0.0000000  0.0000000
#> 12  0.00000000  0.0000000    0.20432063 -0.21789893     0.0000000  0.0000000
#> 13  0.00000000  0.0000000    3.26491288 -1.55732132     0.0000000  0.0000000
#> 14  0.00000000  0.0000000   -6.82284001 -1.10615534     0.0000000  0.0000000
#> 15  0.00000000  0.0000000    0.76837455  1.85727317     0.0000000  0.0000000
#> 16  0.00000000  0.0000000   -7.88216780 10.71574631     0.0000000  0.0000000
#> 17  0.00000000  0.0000000    1.00562559  0.39004188     0.0000000  0.0000000
#> 18  0.00000000  0.0000000   -1.13610174 -0.61873650     0.0000000  0.0000000
#> 19  0.00000000  0.0000000    7.02778962 -7.40123388     0.0000000  0.0000000
#> 20  0.00000000  0.0000000    3.55620531 -2.10187319     0.0000000  0.0000000
#> 21  0.00000000  0.0000000    0.00000000  0.00000000     2.2664960  6.5569956
#> 22  0.00000000  0.0000000    0.00000000  0.00000000     1.5600661 -1.6637416
#> 23  0.00000000  0.0000000    0.00000000  0.00000000    -0.3647111  0.1739625
#> 24  0.00000000  0.0000000    0.00000000  0.00000000    -4.6683977 -0.7568656
#> 25  0.00000000  0.0000000    0.00000000  0.00000000    -1.4640503 -3.5388228
#> 26  0.00000000  0.0000000    0.00000000  0.00000000     5.0226778 -6.8282917
#> 27  0.00000000  0.0000000    0.00000000  0.00000000     2.4031674  0.9320923
#> 28  0.00000000  0.0000000    0.00000000  0.00000000     0.4382241  0.2386628
#> 29  0.00000000  0.0000000    0.00000000  0.00000000    -3.9308824  4.1397625
#> 30  0.00000000  0.0000000    0.00000000  0.00000000    -1.2625898  0.7462459


## Monte Carlo approach, simple trend test example

m <- categorical(lvm(),~x,K=5)
regression(m,additive=TRUE) <- y~x
d <- simulate(m,100,seed=1,'y~x'=0.1)
l <- lm(y~-1+factor(x),data=d)

f <- function(x) coef(lm(x~seq_along(x)))[2]
null <- rep(mean(coef(l)),length(coef(l)))
##  just need to make sure we simulate under H0: slope=0
estimate(l,f,R=1e2,null.sim=null)
#> 100 replications
#> 
#>          seq_along(x)
#> Mean         0.014615
#> SD           0.064600
#>                      
#> 2.5%        -0.094477
#> 97.5%        0.146613
#>                      
#> Estimate     0.080949
#> P-value      0.180000
#> 

estimate(l,f)
#>              Estimate Std.Err     2.5%  97.5% P-value
#> seq_along(x)  0.08095 0.06135 -0.03929 0.2012   0.187

# ------ influence function calculus -------
a <- estimate(coef = c("a" = 0.5), IC = rnorm(10), id = 1:10)
b <- estimate(coef = c("b" = 0.8), IC = rnorm(10), id = 1:10)

e <- c(a, b) # merge
merge(a, b)
#>   Estimate Std.Err    2.5% 97.5% P-value
#> a      0.5  0.2778 -0.0444 1.044 0.07184
#> b      0.8  0.3420  0.1296 1.470 0.01934
c(e1=a, b) # naming of par
#>    Estimate Std.Err    2.5% 97.5% P-value
#> e1      0.5  0.2778 -0.0444 1.044 0.07184
#> b       0.8  0.3420  0.1296 1.470 0.01934
labels(e, c("p1", "p2")) # renaming parameters
#>    Estimate Std.Err    2.5% 97.5% P-value
#> p1      0.5  0.2778 -0.0444 1.044 0.07184
#> p2      0.8  0.3420  0.1296 1.470 0.01934
e["a"] # subset
#>   Estimate Std.Err    2.5% 97.5% P-value
#> a      0.5  0.2778 -0.0444 1.044 0.07184
subset(e, "a")
#>   Estimate Std.Err    2.5% 97.5% P-value
#> a      0.5  0.2778 -0.0444 1.044 0.07184

# pipes
# c(a, b) |>
#  transform(function(x) x^2) |>
#  subset("a") |>
#  labels("sq")

# Parameter transformation with automatic calculation of derivatives
a * b
#>   Estimate Std.Err    2.5%  97.5% P-value
#> a      0.4  0.2799 -0.1486 0.9486   0.153
(3 * cos(a) / sqrt(b) + 1) / a
#>   Estimate Std.Err   2.5% 97.5% P-value
#> a    7.887   5.418 -2.733 18.51  0.1455
expit(c(a,b))
#>   Estimate Std.Err   2.5%  97.5%   P-value
#> a   0.6225 0.06527 0.4945 0.7504 1.484e-21
#> b   0.6900 0.07317 0.5466 0.8334 4.099e-21
c(sum=sum(e), sum2=a+b,
  prod=prod(e), prod2=a*b)
#>       Estimate Std.Err    2.5%  97.5%  P-value
#> sum        1.3  0.4399  0.4379 2.1621 0.003121
#> sum2       1.3  0.4399  0.4379 2.1621 0.003121
#> prod       0.4  0.2799 -0.1486 0.9486 0.153010
#> prod2      0.4  0.2799 -0.1486 0.9486 0.153010
e %*% e # inner prod.
#>    Estimate Std.Err    2.5% 97.5% P-value
#> p1     0.89  0.6128 -0.3112 2.091  0.1464
c(1, 2) %*% e
#>    Estimate Std.Err   2.5% 97.5%  P-value
#> p1      2.1  0.7374 0.6547 3.545 0.004403
c(pow = a^b)
#>     Estimate Std.Err     2.5% 97.5% P-value
#> pow   0.5743  0.2897 0.006493 1.142 0.04744
a^c(0.5, 2)
#>    Estimate Std.Err    2.5%  97.5%   P-value
#> p1   0.7071  0.1964  0.3222 1.0921 0.0003179
#> p2   0.2500  0.2778 -0.2944 0.7944 0.3680875
c(b=e["a"] * e["b"] / a, also.b=e["b"])
#>        Estimate Std.Err   2.5% 97.5% P-value
#> b           0.8   0.342 0.1296  1.47 0.01934
#> also.b      0.8   0.342 0.1296  1.47 0.01934

B <- rbind(c(1,-1), c(1,0), c(0,1))
B %*% e
#>           Estimate Std.Err    2.5%  97.5% P-value
#> [a] - [b]     -0.3  0.4414 -1.1651 0.5651 0.49671
#> a              0.5  0.2778 -0.0444 1.0444 0.07184
#> b              0.8  0.3420  0.1296 1.4704 0.01934
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [a] - [b] = 0
#>   [a] = 0
#>   [b] = 0 
#>  
#> chisq = 8.7407, df = 2, p-value = 0.01265
e == 1 # wald-test, null-hypothesis H0: b=1
#>   Estimate Std.Err    2.5% 97.5% P-value
#> a      0.5  0.2778 -0.0444 1.044 0.07184
#> b      0.8  0.3420  0.1296 1.470 0.55874
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [a] = 1
#>   [b] = 1 
#>  
#> chisq = 3.5899, df = 2, p-value = 0.1661
e == c(1,2)
#>   Estimate Std.Err    2.5% 97.5%  P-value
#> a      0.5  0.2778 -0.0444 1.044 0.071841
#> b      0.8  0.3420  0.1296 1.470 0.000451
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [a] = 1
#>   [b] = 2 
#>  
#> chisq = 15.5935, df = 2, p-value = 0.0004111
B %*% e == 1
#>           Estimate Std.Err    2.5%  97.5%  P-value
#> [a] - [b]     -0.3  0.4414 -1.1651 0.5651 0.003227
#> a              0.5  0.2778 -0.0444 1.0444 0.071841
#> b              0.8  0.3420  0.1296 1.4704 0.558742
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[a] - [b]] = 1
#>   [a] = 1
#>   [b] = 1 
#>  
#> chisq = 9.145, df = 2, p-value = 0.01033
```
