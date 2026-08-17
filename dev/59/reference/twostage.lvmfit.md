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
m <- lvm(c(x1,x2,x3)~f1,f1~z,
         c(y1,y2,y3)~f2,f2~f1+z)
latent(m) <- ~f1+f2
d <- simulate(m,100,p=c("f2,f2"=2,"f1,f1"=0.5),seed=1)

## Full MLE
ee <- estimate(m,d)

## Manual two-stage
if (FALSE) { # \dontrun{
m1 <- lvm(c(x1,x2,x3)~f1,f1~z); latent(m1) <- ~f1
e1 <- estimate(m1,d)
pp1 <- predict(e1,f1~x1+x2+x3)

d$u1 <- pp1[,]
d$u2 <- pp1[,]^2+attr(pp1,"cond.var")[1]
m2 <- lvm(c(y1,y2,y3)~eta,c(y1,eta)~u1+u2+z); latent(m2) <- ~eta
e2 <- estimate(m2,d)
} # }

## Two-stage
m1 <- lvm(c(x1,x2,x3)~f1,f1~z); latent(m1) <- ~f1
m2 <- lvm(c(y1,y2,y3)~eta,c(y1,eta)~u1+u2+z); latent(m2) <- ~eta
pred <- function(mu,var,data,...)
    cbind("u1"=mu[,1],"u2"=mu[,1]^2+var[1])
(mm <- twostage(m1,model2=m2,data=d,predict.fun=pred))
#>                     Estimate Std. Error  Z-value   P-value
#> Measurements:                                             
#>    y2~eta            0.96270    0.12462  7.72525    <1e-12
#>    y3~eta            0.97886    0.12477  7.84516    <1e-12
#> Regressions:                                              
#>    y1~u1            -0.10923    0.24754 -0.44125     0.659
#>    y1~u2            -0.00916    0.02941 -0.31149    0.7554
#>    y1~z             -0.09088    0.25697 -0.35364    0.7236
#>     eta~u1           1.23462    0.24891  4.96011 7.045e-07
#>     eta~u2           0.00912    0.02417  0.37737    0.7059
#>     eta~z            0.84531    0.27943  3.02508  0.002486
#> Intercepts:                                               
#>    y2               -0.19048    0.16526 -1.15257    0.2491
#>    y3                0.00979    0.18282  0.05354    0.9573
#>    eta              -0.16805    0.22950 -0.73226     0.464
#> Residual Variances:                                       
#>    y1                1.15103    0.24219  4.75250          
#>    y2                0.97707    0.20212  4.83411          
#>    y3                1.13661    0.20506  5.54289          
#>    eta               1.58985    0.37736  4.21312          

if (interactive()) {
    pf <- function(p) p["eta"]+p["eta~u1"]*u + p["eta~u2"]*u^2
    plot(mm,f=pf,data=data.frame(u=seq(-2,2,length.out=100)),lwd=2)
}

 ## Reduce test timing
## Splines
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

kn <- c(-1,0,1)
nonlinear(m2,type="spline",knots=kn) <- eta2~eta1
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
 ## Reduce test timing

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
#>                     Estimate Std. Error  Z-value  P-value
#> Measurements:                                            
#>    y2~eta2           1.03076    0.03876 26.59191   <1e-12
#>    y3~eta2           0.99635    0.04036 24.68596   <1e-12
#> Regressions:                                             
#>    eta2~eta1_1       2.47344    0.22699 10.89681   <1e-12
#>    eta2~eta1_2      -0.48653    0.05555 -8.75771   <1e-12
#> Intercepts:                                              
#>    y2                0.02005    0.06568  0.30528   0.7602
#>    y3                0.07986    0.06779  1.17797   0.2388
#>    eta2              1.16241    0.16093  7.22290   <1e-12
#> Residual Variances:                                      
#>    y1                1.11229    0.10786 10.31201         
#>    y2                1.00927    0.09818 10.28011         
#>    y3                1.09169    0.09757 11.18840         
#>    eta2              1.81941    0.18188 10.00346         
 # \dontrun{}
```
