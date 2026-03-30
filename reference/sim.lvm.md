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
#>   y          x z
#> 1 0 -1.0871117 0
#> 2 0 -0.5961507 0
#> 3 1  1.2456606 1
#> 4 0 -1.9011011 0
#> 5 1  0.9317311 0
#> 6 0 -1.0195031 0

e <- estimate(m,d,estimator="glm")
e
#>              Estimate Std. Error  Z-value   P-value
#> Regressions:                                       
#>    y~x        0.97832    0.09345 10.46919    <1e-12
#>    y~z        1.04098    0.17281  6.02404 1.701e-09
#>     x~z       1.11747    0.06529 17.11593    <1e-12
#> Intercepts:                                        
#>    y         -0.01437    0.09874 -0.14555    0.8843
#>    x         -0.05580    0.04649 -1.20022    0.2301
#> Dispersion:                                        
#>    x          1.06783                              
## Simulate a few observation from estimated model
sim(e,n=5)
#>   y          x z
#> 1 0 -0.1463386 0
#> 2 1  0.9463616 1
#> 3 1  0.1632782 1
#> 4 0 -1.2455533 0
#> 5 1  2.8957783 1

##################################################
## Poisson
##################################################
distribution(m,~y) <- poisson.lvm()
d <- sim(m,1e4,p=c(y=-1,"y~x"=2,z=1))
head(d)
#>    y          x z
#> 1  1 -0.1831437 1
#> 2  3  0.4235717 1
#> 3 19  1.4638566 1
#> 4  0 -1.3786736 1
#> 5  0  0.7720821 0
#> 6  2 -0.4325138 1
estimate(m,d,estimator="glm")
#>                Estimate Std. Error    Z-value  P-value
#> Regressions:                                          
#>    y~x          2.00013    0.00170 1179.77248   <1e-12
#>    y~z          0.99804    0.01139   87.63590   <1e-12
#>     x~z         0.99225    0.02256   43.99038   <1e-12
#> Intercepts:                                           
#>    y           -0.99834    0.01161  -85.98726   <1e-12
#>    x            0.00295    0.01927    0.15321   0.8782
#> Dispersion:                                           
#>    x            1.00388                               
mean(d$z); lava:::expit(1)
#> [1] 0.7255
#> [1] 0.7310586

summary(lm(y~x,sim(lvm(y[1:2]~4*x),1e3)))
#> 
#> Call:
#> lm(formula = y ~ x, data = sim(lvm(y[1:2] ~ 4 * x), 1000))
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -5.3176 -0.9817  0.0425  0.9454  5.3563 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  1.01691    0.04476   22.72   <2e-16 ***
#> x            3.98601    0.04489   88.80   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 1.414 on 998 degrees of freedom
#> Multiple R-squared:  0.8877, Adjusted R-squared:  0.8875 
#> F-statistic:  7885 on 1 and 998 DF,  p-value: < 2.2e-16
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
#> (Intercept) 0.503872   0.005074   99.30   <2e-16 ***
#> x           1.015755   0.016366   62.07   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for Gamma family taken to be 0.5154657)
#> 
#>     Null deviance: 8439.4  on 9999  degrees of freedom
#> Residual deviance: 5553.3  on 9998  degrees of freedom
#> AIC: 20719
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
#>  [1] -0.39232531 -0.65344132 -0.25008661  1.02265025 -0.29470894 -0.75875908
#>  [7] -0.61782561 -0.01173492  0.22658620 -0.85462986

##################################################
### Beta
##################################################
m <- lvm()
distribution(m,~y) <- beta.lvm(alpha=2,beta=1)
var(sim(m,100,"y,y"=2))
#>          y
#> y 1.038929
distribution(m,~y) <- beta.lvm(alpha=2,beta=1,scale=FALSE)
var(sim(m,100))
#>            y
#> y 0.06458889

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
#> -3.1733 -0.6659 -0.0400  0.7135  3.6610 
#> 
#> Coefficients:
#>                Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)     0.03640    0.06153   0.592    0.554    
#> x               0.97845    0.04709  20.778   <2e-16 ***
#> z               1.01007    0.05412  18.663   <2e-16 ***
#> I(z > 0)TRUE   -0.01420    0.10558  -0.134    0.893    
#> x:I(z > 0)TRUE  1.02565    0.06515  15.744   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 1.015 on 995 degrees of freedom
#> Multiple R-squared:  0.7725, Adjusted R-squared:  0.7716 
#> F-statistic: 844.9 on 4 and 995 DF,  p-value: < 2.2e-16
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
#>             y            x          z           t
#> 1   0.5579720  2.120772086 0.33831431  0.18876989
#> 2   1.4494352 -0.135160928 3.33297133  4.83092594
#> 3   0.4875242  1.070643367 0.04153628  0.02024994
#> 4   0.2123152  1.886413966 2.46505719  0.52336918
#> 5   0.2478579  1.136382241 0.34452875  0.08539416
#> 6  -0.7640068 -1.826504666 0.01332605 -0.01018119
#> 7   0.5955366  1.284139328 0.23062140  0.13734349
#> 8   0.7351995 -0.073520000 0.44168647  0.32472769
#> 9   0.9358395 -0.001793702 0.57737245  0.54032796
#> 10  0.9575549  0.155946016 1.98820808  1.90381833
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
#>          y1       y2       y3
#> y1 2.673436 1.237220 1.280251
#> y2 1.237220 1.979091 1.121231
#> y3 1.280251 1.121231 2.231319

## Syntax also useful for univariate generators, e.g.
m <- lvm(y~x+z)
distribution(m,~y,TRUE) <- function(n) rnorm(n,mean=1000)
sim(m,5)
#>           y          x           z
#> 1  999.6598 -1.2941958 -0.21599390
#> 2  995.4996 -1.3899030 -1.79740633
#> 3 1002.6299  1.0073420  0.46536012
#> 4 1000.0435 -1.4961390  0.72058077
#> 5 1000.0448  0.5213245 -0.00315081
distribution(m,~y,"m1",0) <- rnorm
sim(m,5)
#>             y          x           z
#> 1 -0.27522454 -1.2971259  0.21634297
#> 2 -1.80708496  0.4350938 -1.58337733
#> 3  0.05877631  0.5108942 -0.38597289
#> 4  2.84314360  1.5773881 -0.02885458
#> 5 -0.35209704 -1.0659894  1.81281571
sim(m,5,p=c(m1=100))
#>          y           x          z
#> 1 101.0689 -0.35860612  0.2469560
#> 2 101.9527  1.93296840 -1.6293429
#> 3 100.7706  0.10528753 -0.4887392
#> 4 100.4315 -0.09999915  1.4974758
#> 5 100.7946 -0.13847376 -0.2257985

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
#> 5006 3017 1977 
glm(y~v, data=sim(m,1e4))
#> 
#> Call:  glm(formula = y ~ v, data = sim(m, 10000))
#> 
#> Coefficients:
#> (Intercept)           vB           vC  
#>     0.02251      2.00064     -1.00417  
#> 
#> Degrees of Freedom: 9999 Total (i.e. Null);  9997 Residual
#> Null Deviance:       22620 
#> Residual Deviance: 9990  AIC: 28380
glm(y~v, data=sim(m,1e4,p=c("y~v:1"=3)))
#> 
#> Call:  glm(formula = y ~ v, data = sim(m, 10000, p = c(`y~v:1` = 3)))
#> 
#> Coefficients:
#> (Intercept)           vB           vC  
#>    0.003158     2.986700    -1.011000  
#> 
#> Degrees of Freedom: 9999 Total (i.e. Null);  9997 Residual
#> Null Deviance:       34540 
#> Residual Deviance: 10350     AIC: 28730

transform(m,v2~v) <- function(x) x=='A'
sim(m,10)
#>    v           y          z    v2
#> 1  B  2.46056403  2.8858782 FALSE
#> 2  B  2.13354634  4.1726117 FALSE
#> 3  B  2.83047801  5.0140621 FALSE
#> 4  A  0.61279199  0.1483464  TRUE
#> 5  A -0.10165782  0.5226807  TRUE
#> 6  A  1.00767457 -0.6966552  TRUE
#> 7  A -0.66799846 -0.9295116  TRUE
#> 8  A  0.87462526 -1.6325467  TRUE
#> 9  A  0.00129238 -1.1372980  TRUE
#> 10 C -2.20096227  0.5693128 FALSE

##################################################
### Pre-calculate object
##################################################
m <- lvm(y~x)
m2 <- sim(m,'y~x'=2)
sim(m,10,'y~x'=2)
#>             y          x
#> 1   1.2233025  0.5982027
#> 2  -1.9169139 -0.1868558
#> 3   1.6462627 -0.1377385
#> 4  -1.2117513 -0.4254388
#> 5   3.0226593  1.6894790
#> 6  -3.9541658 -1.5825774
#> 7   0.8057423  0.2298156
#> 8  -2.1578435 -1.3413077
#> 9   1.8430087  0.6336803
#> 10  1.1745427 -0.2815957
sim(m2,10) ## Faster
#>            y          x
#> 1  -4.625097 -2.2077901
#> 2  -3.239577 -2.4106877
#> 3   2.582984  0.8197093
#> 4   2.435479  0.9182792
#> 5   3.031488  1.2826166
#> 6  -2.039096 -1.3749699
#> 7  -3.497009 -1.8324686
#> 8   2.049103  1.1580864
#> 9  -1.319555  0.2221127
#> 10 -5.157303 -1.4383467
```
