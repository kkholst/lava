# Simulate model

Simulate data from a general SEM model including non-linear effects and
general link and distribution of variables.

## Usage

``` r
# S3 method for class 'lvm'
sim(x, n = NULL, p = NULL, normal = FALSE, cond = FALSE,
sigma = 1, rho = 0.5, X = NULL, unlink=FALSE, latent=TRUE,
use.labels = TRUE, seed=NULL, ...)
```

## Arguments

- x:

  Model object

- n:

  Number of simulated values/individuals

- p:

  Parameter value (optional)

- normal:

  Logical indicating whether to simulate data from a multivariate normal
  distribution conditional on exogenous variables hence ignoring
  functional/distribution definition

- cond:

  for internal use

- sigma:

  Default residual variance (1)

- rho:

  Default covariance parameter (0.5)

- X:

  Optional matrix of fixed values of variables (manipulation)

- unlink:

  Return Inverse link transformed data

- latent:

  Include latent variables (default TRUE)

- use.labels:

  convert categorical variables to factors before applying
  transformation

- seed:

  Random seed

- ...:

  Additional arguments to be passed to the low level functions

## Author

Klaus K. Holst

## Examples

``` r
##################################################
## Logistic regression
##################################################
m <- lvm(y~x+z)
regression(m) <- x~z
distribution(m,~y+z) <- binomial.lvm("logit")
d <- sim(m,1e3)
head(d)
#>   y           x z
#> 1 1  0.40384929 1
#> 2 0  0.24566058 0
#> 3 0 -1.90110113 0
#> 4 1  0.93173113 0
#> 5 1 -0.01950306 1
#> 6 0 -0.39374796 0

e <- estimate(m,d,estimator="glm")
e
#>              Estimate Std. Error  Z-value   P-value
#> Regressions:                                       
#>    y~x        0.97754    0.09021 10.83680    <1e-12
#>    y~z        1.10148    0.17286  6.37194 1.866e-10
#>     x~z       0.94600    0.06534 14.47800    <1e-12
#> Intercepts:                                        
#>    y          0.11229    0.09813  1.14427    0.2525
#>    x          0.03084    0.04550  0.67781    0.4979
#> Dispersion:                                        
#>    x          1.06937                              
## Simulate a few observation from estimated model
sim(e,n=5)
#>   y          x z
#> 1 1  0.8614565 1
#> 2 0 -0.8681925 0
#> 3 0 -0.2137645 1
#> 4 1  2.8122735 1
#> 5 1  1.4170839 1

##################################################
## Poisson
##################################################
distribution(m,~y) <- poisson.lvm()
d <- sim(m,1e4,p=c(y=-1,"y~x"=2,z=1))
head(d)
#>    y           x z
#> 1  2  0.42357168 1
#> 2 15  1.46385660 1
#> 3  0 -2.37867359 0
#> 4 45  1.77208207 1
#> 5  0 -0.43251382 1
#> 6  0 -0.02905308 0
estimate(m,d,estimator="glm")
#>                Estimate Std. Error    Z-value  P-value
#> Regressions:                                          
#>    y~x          2.00130    0.00168 1188.26216   <1e-12
#>    y~z          0.98389    0.01204   81.73191   <1e-12
#>     x~z         1.03013    0.02254   45.70045   <1e-12
#> Intercepts:                                           
#>    y           -0.98389    0.01220  -80.66118   <1e-12
#>    x           -0.02456    0.01925   -1.27584    0.202
#> Dispersion:                                           
#>    x            1.00378                               
mean(d$z); lava:::expit(1)
#> [1] 0.7254
#> [1] 0.7310586

summary(lm(y~x,sim(lvm(y[1:2]~4*x),1e3)))
#> 
#> Call:
#> lm(formula = y ~ x, data = sim(lvm(y[1:2] ~ 4 * x), 1000))
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -4.2200 -1.0532  0.0316  0.9507  4.5948 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.99506    0.04544   21.90   <2e-16 ***
#> x            4.00185    0.04426   90.41   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 1.437 on 998 degrees of freedom
#> Multiple R-squared:  0.8912, Adjusted R-squared:  0.8911 
#> F-statistic:  8175 on 1 and 998 DF,  p-value: < 2.2e-16
#> 

##################################################
### Gamma distribution
##################################################
m <- lvm(y~x)
distribution(m,~y+x) <- list(Gamma.lvm(shape=2),binomial.lvm())
intercept(m,~y) <- 0.5
d <- sim(m,1e4)
summary(g <- glm(y~x,family=Gamma(),data=d))
#> 
#> Call:
#> glm(formula = y ~ x, family = Gamma(), data = d)
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) 0.507575   0.005117   99.20   <2e-16 ***
#> x           1.002598   0.016261   61.66   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for Gamma family taken to be 0.5151108)
#> 
#>     Null deviance: 8368.5  on 9999  degrees of freedom
#> Residual deviance: 5546.7  on 9998  degrees of freedom
#> AIC: 20669
#> 
#> Number of Fisher Scoring iterations: 6
#> 
if (FALSE) MASS::gamma.shape(g) # \dontrun{}

args(lava::Gamma.lvm)
#> function (link = "inverse", shape, rate, unit = FALSE, var = FALSE, 
#>     log = FALSE, ...) 
#> NULL
distribution(m,~y) <- Gamma.lvm(shape=2,log=TRUE)
sim(m,10,p=c(y=0.5))[,"y"]
#>  [1] -1.3148925  0.0231124  0.6164099 -0.3583493 -0.0718344 -0.2697476
#>  [7] -1.1070485  1.0313895  0.7453192 -1.1577031

##################################################
### Beta
##################################################
m <- lvm()
distribution(m,~y) <- beta.lvm(alpha=2,beta=1)
var(sim(m,100,"y,y"=2))
#>          y
#> y 1.027512
distribution(m,~y) <- beta.lvm(alpha=2,beta=1,scale=FALSE)
var(sim(m,100))
#>            y
#> y 0.04661863

##################################################
### Transform
##################################################
m <- lvm()
transform(m,xz~x+z) <- function(x) x[1]*(x[2]>0)
regression(m) <- y~x+z+xz
d <- sim(m,1e3)
summary(lm(y~x+z + x*I(z>0),d))
#> 
#> Call:
#> lm(formula = y ~ x + z + x * I(z > 0), data = d)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -3.6128 -0.6418  0.0358  0.6943  3.0594 
#> 
#> Coefficients:
#>                  Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)    -0.0103490  0.0627031  -0.165    0.869    
#> x               1.0236548  0.0452818  22.606   <2e-16 ***
#> z               1.0021651  0.0543581  18.436   <2e-16 ***
#> I(z > 0)TRUE    0.0001651  0.1087734   0.002    0.999    
#> x:I(z > 0)TRUE  0.9391138  0.0626959  14.979   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 1.011 on 995 degrees of freedom
#> Multiple R-squared:  0.7824, Adjusted R-squared:  0.7815 
#> F-statistic: 894.4 on 4 and 995 DF,  p-value: < 2.2e-16
#> 

##################################################
### Non-random variables
##################################################
m <- lvm()
distribution(m,~x+z+v+w) <- list(Sequence.lvm(0,5),## Seq. 0 to 5 by 1/n
                               Binary.lvm(),       ## Vector of ones
                               Binary.lvm(0.5),    ##  0.5n 0, 0.5n 1
                               Binary.lvm(interval=list(c(0.3,0.5),c(0.8,1))))
sim(m,10)
#>            x z v w
#> 1  0.0000000 1 0 0
#> 2  0.5555556 1 0 0
#> 3  1.1111111 1 0 1
#> 4  1.6666667 1 0 1
#> 5  2.2222222 1 0 1
#> 6  2.7777778 1 1 0
#> 7  3.3333333 1 1 0
#> 8  3.8888889 1 1 1
#> 9  4.4444444 1 1 1
#> 10 5.0000000 1 1 1

##################################################
### Cox model
### piecewise constant hazard
################################################
m <- lvm(t~x)
rates <- c(1,0.5); cuts <- c(0,5)
## Constant rate: 1 in [0,5), 0.5 in [5,Inf)
distribution(m,~t) <- coxExponential.lvm(rate=rates,timecut=cuts)

if (FALSE) { # \dontrun{
    d <- sim(m,2e4,p=c("t~x"=0.1)); d$status <- TRUE
    plot(timereg::aalen(survival::Surv(t,status)~x,data=d,
                        resample.iid=0,robust=0),spec=1)
    L <- approxfun(c(cuts,max(d$t)),f=1,
                   cumsum(c(0,rates*diff(c(cuts,max(d$t))))),
                   method="linear")
    curve(L,0,100,add=TRUE,col="blue")
} # }

##################################################
### Cox model
### piecewise constant hazard, gamma frailty
##################################################
m <- lvm(y~x+z)
rates <- c(0.3,0.5); cuts <- c(0,5)
distribution(m,~y+z) <- list(coxExponential.lvm(rate=rates,timecut=cuts),
                             loggamma.lvm(rate=1,shape=1))
if (FALSE) { # \dontrun{
    d <- sim(m,2e4,p=c("y~x"=0,"y~z"=0)); d$status <- TRUE
    plot(timereg::aalen(survival::Surv(y,status)~x,data=d,
                        resample.iid=0,robust=0),spec=1)
    L <- approxfun(c(cuts,max(d$y)),f=1,
                   cumsum(c(0,rates*diff(c(cuts,max(d$y))))),
                   method="linear")
    curve(L,0,100,add=TRUE,col="blue")
} # }
## Equivalent via transform (here with Aalens additive hazard model)
m <- lvm(y~x)
distribution(m,~y) <- aalenExponential.lvm(rate=rates,timecut=cuts)
distribution(m,~z) <- Gamma.lvm(rate=1,shape=1)
transform(m,t~y+z) <- prod
sim(m,10)
#>             y          x           z            t
#> 1  -3.2840809 -0.5313786 0.081230892 -0.266768824
#> 2   0.4684835  1.0017989 0.842582853  0.394736176
#> 3   0.1286571  0.7976910 1.952164664  0.251159881
#> 4   3.8958800  0.3282328 2.346598408  9.142065876
#> 5  -0.9036212 -0.4915184 0.087873678 -0.079404519
#> 6   0.1126743  1.0220275 1.532014757  0.172618715
#> 7   0.6427252  0.3751732 0.508685001  0.326944669
#> 8   1.4215918 -0.1591446 0.590513776  0.839469548
#> 9   0.3313485  1.3396558 0.205874305  0.068216143
#> 10  0.3355325  0.1197760 0.009979276  0.003348371
## Shared frailty
m <- lvm(c(t1,t2)~x+z)
rates <- c(1,0.5); cuts <- c(0,5)
distribution(m,~y) <- aalenExponential.lvm(rate=rates,timecut=cuts)
distribution(m,~z) <- loggamma.lvm(rate=1,shape=1)
if (FALSE) { # \dontrun{
mets::fast.reshape(sim(m,100),varying="t")
} # }

##################################################
### General multivariate distributions
##################################################
if (FALSE) { # \dontrun{
m <- lvm()
distribution(m,~y1+y2,oratio=4) <- VGAM::rbiplackcop
ksmooth2(sim(m,1e4),rgl=FALSE,theta=-20,phi=25)

m <- lvm()
distribution(m,~z1+z2,"or1") <- VGAM::rbiplackcop
distribution(m,~y1+y2,"or2") <- VGAM::rbiplackcop
sim(m,10,p=c(or1=0.1,or2=4))
} # }

m <- lvm()
distribution(m,~y1+y2+y3,TRUE) <- function(n,...) rmvn0(n,sigma=diag(3)+1)
var(sim(m,100))
#>           y1        y2        y3
#> y1 1.7249852 0.6674309 0.5390163
#> y2 0.6674309 1.4423047 0.7318100
#> y3 0.5390163 0.7318100 1.6936224

## Syntax also useful for univariate generators, e.g.
m <- lvm(y~x+z)
distribution(m,~y,TRUE) <- function(n) rnorm(n,mean=1000)
sim(m,5)
#>           y          x          z
#> 1  994.5267 -2.3053292 -1.8824933
#> 2  998.4568 -1.3360425 -0.3062621
#> 3 1002.8650 -0.3013739  1.0649400
#> 4 1001.3514  1.2848582  0.9951471
#> 5 1000.6749  1.5686740 -0.4457475
distribution(m,~y,"m1",0) <- rnorm
sim(m,5)
#>            y          x          z
#> 1 -0.3940713 -0.6209352 -0.6273853
#> 2 -0.7275340  0.5040422  0.3378076
#> 3 -1.3146335 -1.6480163  0.7666236
#> 4 -1.9969884 -0.9574870 -1.0264748
#> 5  2.5320341  2.2250715  0.4485908
sim(m,5,p=c(m1=100))
#>           y          x          z
#> 1 101.46870  1.6884779 -0.5837411
#> 2  95.68423 -1.1575305 -1.1663874
#> 3 103.15080  0.5007453  1.5080476
#> 4  97.79707 -1.8126868 -0.9694514
#> 5  95.86543 -1.4649504 -0.9732468

##################################################
### Regression design in other parameters
##################################################
## Variance heterogeneity
m <- lvm(y~x)
distribution(m,~y) <- function(n,mean,x) rnorm(n,mean,exp(x)^.5)
if (interactive()) plot(y~x,sim(m,1e3))
## Alternaively, calculate the standard error directly
addvar(m) <- ~sd ## If 'sd' should be part of the resulting data.frame
constrain(m,sd~x) <- function(x) exp(x)^.5
distribution(m,~y) <- function(n,mean,sd) rnorm(n,mean,sd)
if (interactive()) plot(y~x,sim(m,1e3))

## Regression on variance parameter
m <- lvm()
regression(m) <- y~x
regression(m) <- v~x
##distribution(m,~v) <- 0 # No stochastic term
## Alternative:
## regression(m) <- v[NA:0]~x
distribution(m,~y) <- function(n,mean,v) rnorm(n,mean,exp(v)^.5)
if (interactive()) plot(y~x,sim(m,1e3))

## Regression on shape parameter in Weibull model
m <- lvm()
regression(m) <- y ~ z+v
regression(m) <- s ~ exp(0.6*x-0.5*z)
distribution(m,~x+z) <- binomial.lvm()
distribution(m,~cens) <- coxWeibull.lvm(scale=1)
distribution(m,~y) <- coxWeibull.lvm(scale=0.1,shape=~s)
eventTime(m) <- time ~ min(y=1,cens=0)

if (interactive()) {
    d <- sim(m,1e3)
    require(survival)
    (cc <- coxph(Surv(time,status)~v+strata(x,z),data=d))
    plot(survfit(cc) ,col=1:4,mark.time=FALSE)
}

##################################################
### Categorical predictor
##################################################
m <- lvm()
## categorical(m,K=3) <- "v"
categorical(m,labels=c("A","B","C")) <- "v"

regression(m,additive=FALSE) <- y~v
if (FALSE) { # \dontrun{
plot(y~v,sim(m,1000,p=c("y~v:2"=3)))
} # }

m <- lvm()
categorical(m,labels=c("A","B","C"),p=c(0.5,0.3)) <- "v"
regression(m,additive=FALSE,beta=c(0,2,-1)) <- y~v
## equivalent to:
## regression(m,y~v,additive=FALSE) <- c(0,2,-1)
regression(m,additive=FALSE,beta=c(0,4,-1)) <- z~v
table(sim(m,1e4)$v)
#> 
#>    A    B    C 
#> 5001 3050 1949 
glm(y~v, data=sim(m,1e4))
#> 
#> Call:  glm(formula = y ~ v, data = sim(m, 10000))
#> 
#> Coefficients:
#> (Intercept)           vB           vC  
#>   -0.003607     2.032486    -0.983431  
#> 
#> Degrees of Freedom: 9999 Total (i.e. Null);  9997 Residual
#> Null Deviance:       22410 
#> Residual Deviance: 9816  AIC: 28200
glm(y~v, data=sim(m,1e4,p=c("y~v:1"=3)))
#> 
#> Call:  glm(formula = y ~ v, data = sim(m, 10000, p = c(`y~v:1` = 3)))
#> 
#> Coefficients:
#> (Intercept)           vB           vC  
#>   -0.001448     3.033790    -0.987477  
#> 
#> Degrees of Freedom: 9999 Total (i.e. Null);  9997 Residual
#> Null Deviance:       34140 
#> Residual Deviance: 9790  AIC: 28170

transform(m,v2~v) <- function(x) x=='A'
sim(m,10)
#>    v          y          z    v2
#> 1  C -1.4228910 -2.3743813 FALSE
#> 2  B  3.7990113  2.3934600 FALSE
#> 3  A  0.7750815 -0.3692173  TRUE
#> 4  C  0.1102358 -0.8214093 FALSE
#> 5  B  1.5153310  6.1654183 FALSE
#> 6  B  2.8497096  5.6700463 FALSE
#> 7  B  2.0529501  4.3634553 FALSE
#> 8  B  2.3471577  3.7715189 FALSE
#> 9  C -2.3345083 -2.2933131 FALSE
#> 10 A  2.1200748 -0.2894096  TRUE

##################################################
### Pre-calculate object
##################################################
m <- lvm(y~x)
m2 <- sim(m,'y~x'=2)
sim(m,10,'y~x'=2)
#>             y          x
#> 1   4.2044588  2.1870327
#> 2  -1.2835247 -0.3542078
#> 3  -0.6652916  0.1662725
#> 4   4.3367737  0.5560635
#> 5  -2.3096870 -1.3264723
#> 6   1.0442205  0.1227436
#> 7  -2.4180767 -0.5474568
#> 8   0.8563611  1.1660111
#> 9  -0.8626151 -0.9314213
#> 10  4.0543733  1.3070954
sim(m2,10) ## Faster
#>             y           x
#> 1   1.3121923  1.50961388
#> 2   1.0176849  0.14595064
#> 3  -3.8403860 -1.48995665
#> 4   0.9179312 -0.13233502
#> 5   0.8011242  0.01558654
#> 6   2.1012753  1.01420451
#> 7  -1.4526002 -1.29120883
#> 8   2.1422520  1.23362195
#> 9   2.8050867  0.67871130
#> 10  2.7790760  1.50309619
```
