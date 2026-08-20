# Two-stage estimator (non-linear SEM)

Two-stage estimator for non-linear structural equation models

## Usage

``` r
# S3 method for class 'lvmfit'
twostage(
  object,
  model2,
  data = NULL,
  predict.fun = NULL,
  id1 = NULL,
  id2 = NULL,
  all = FALSE,
  formula = NULL,
  std.err = TRUE,
  ...
)
```

## Arguments

- object:

  Stage 1 measurement model

- model2:

  Stage 2 SEM

- data:

  data.frame

- predict.fun:

  Prediction of latent variable

- id1:

  Optional id-variable (stage 1 model)

- id2:

  Optional id-variable (stage 2 model)

- all:

  If TRUE return additional output (naive estimates)

- formula:

  optional formula specifying non-linear relation

- std.err:

  If FALSE calculations of standard errors will be skipped

- ...:

  Additional arguments to lower level functions

## Examples

``` r
if (FALSE)  ## Reduce test timing
# simulated example (only linear effects)
m <- lvm(c(x1,x2,x3)~f1,f1~z,
         c(y1,y2,y3)~f2,f2~f1+z)
latent(m) <- ~f1+f2
#> Error: object 'm' not found
d <- simulate(m,1000,p=c("f2,f2"=2,"f1,f1"=0.5),seed=1)
#> Error: object 'm' not found

## Full MLE
ee <- estimate(m,d)
#> Error: object 'm' not found

## Manual two-stage
m1 <- lvm(c(x1,x2,x3)~f1,f1~z); latent(m1) <- ~f1
e1 <- estimate(m1,d)
#> Error: object 'd' not found
pp1 <- predict(e1,f1~x1+x2+x3)
#> Error: object 'e1' not found

d$u1 <- pp1[,]
#> Error: object 'pp1' not found
d$u2 <- pp1[,]^2+attr(pp1,"cond.var")[1]
#> Error: object 'pp1' not found
m2 <- lvm(c(y1,y2,y3)~eta,c(y1,eta)~u1+u2+z); latent(m2) <- ~eta
e2 <- estimate(m2,d)
#> Error: object 'd' not found

## twostage method:
m1 <- lvm(c(x1,x2,x3)~f1,f1~z); latent(m1) <- ~f1
m2 <- lvm(c(y1,y2,y3)~eta,c(y1,eta)~u1+u2+z); latent(m2) <- ~eta
pred <- function(mu,var,data,...)
    cbind("u1"=mu[,1],"u2"=mu[,1]^2+var[1])
(mm <- twostage(m1,model2=m2,data=d,predict.fun=pred))
#> Error: object 'd' not found

  if (interactive()) {
    pf <- function(p) p["eta"]+p["eta~u1"]*u + p["eta~u2"]*u^2
    plot(mm,f=pf,data=data.frame(u=seq(-2,2,length.out=100)),lwd=2)
  }
 # \dontrun{}

## Quadratic example
f <- function(x) cos(2*x)+x+-0.25*x^2
m <- lvm(x1+x2+x3~eta1, y1+y2+y3~eta2, latent=~eta1+eta2)
functional(m, eta2~eta1) <- f
d <- sim(m,500,seed=1,latent=TRUE)
m1 <- lvm(x1+x2+x3~eta1,latent=~eta1)
m2 <- lvm(y1+y2+y3~eta2,latent=~eta2)
mm <- twostage(m1,m2,formula=eta2~eta1,type="spline")
if (interactive()) plot(mm)

nonlinear(m2,type="quadratic") <- eta2~eta1
a <- twostage(m1,m2,data=d)
if (interactive()) plot(a)
if (FALSE)  ## Reduce test timing##'
## Splines
kn <- c(-1,0,1)
nonlinear(m2,type="spline",knots=kn) <- eta2~eta1
#> Error: object 'kn' not found
a <- twostage(m1,m2,data=d)
x <- seq(-3,3,by=0.1)
y <- predict(a, newdata=data.frame(eta1=x))

if (interactive()) {
  plot(eta2~eta1, data=d)
  lines(x,y, col="red", lwd=5)

  p <- estimate(a,f=function(p) predict(a,p=p,newdata=x))$coefmat
  plot(eta2~eta1, data=d)
  lines(x,p[,1], col="red", lwd=5)
  confband(x,lower=p[,3],upper=p[,4],center=p[,1], polygon=TRUE, col=Col(2,0.2))

  l1 <- lm(eta2~splines::ns(eta1,knots=kn),data=d)
  p1 <- predict(l1,newdata=data.frame(eta1=x),interval="confidence")
  lines(x,p1[,1],col="green",lwd=5)
  confband(x,lower=p1[,2],upper=p1[,3],center=p1[,1], polygon=TRUE, col=Col(3,0.2))
}
 # \dontrun{} ## Reduce test timing

if (FALSE)  ## Reduce timing
 ## Cross-validation example
 ma <- lvm(c(x1,x2,x3)~u,latent=~u)
 ms <- functional(ma, y~u, value=function(x) -.4*x^2)
#> Error: object 'ma' not found
 d <- sim(ms,500)#,seed=1)
#> Error: object 'ms' not found
 ea <- estimate(ma,d)
#> Error: object 'ma' not found

 mb <- lvm()
 mb1 <- nonlinear(mb,type="linear",y~u)
 mb2 <- nonlinear(mb,type="quadratic",y~u)
 mb3 <- nonlinear(mb,type="spline",knots=c(-3,-1,0,1,3),y~u)
 mb4 <- nonlinear(mb,type="spline",knots=c(-3,-2,-1,0,1,2,3),y~u)
 ff <- lapply(list(mb1,mb2,mb3,mb4),
      function(m) function(data,...) twostage(ma,m,data=data,st.derr=FALSE))
 a <- cv(ff,data=d,rep=1)
#> Error in cv(ff, data = d, rep = 1): could not find function "cv"
 a
#>                     Estimate Std. Error  Z-value   P-value
#> Measurements:                                             
#>    y2~eta2           1.03159    0.03878 26.60368    <1e-12
#>    y3~eta2           0.99823    0.04055 24.61738    <1e-12
#> Regressions:                                              
#>    eta2~eta1_1       1.02678    0.09758 10.52288    <1e-12
#>    eta2~eta1_2      -0.45396    0.07186 -6.31753 2.658e-10
#> Intercepts:                                               
#>    y2                0.02026    0.06571  0.30828    0.7579
#>    y3                0.08032    0.06785  1.18390    0.2365
#>    eta2              0.34096    0.10913  3.12452  0.001781
#> Residual Variances:                                       
#>    y1                1.11844    0.10884 10.27618          
#>    y2                1.01001    0.09859 10.24506          
#>    y3                1.08519    0.09822 11.04895          
#>    eta2              1.89819    0.18571 10.22130          
 # \dontrun{}
```
