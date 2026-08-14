# Two-stage (non-linear) measurement error

Two-stage measurement error

## Usage

``` r
measurement.error(
  model1,
  formula,
  data = parent.frame(),
  predictfun = function(mu, var, data, ...) mu[, 1]^2 + var[1],
  id1,
  id2,
  ...
)
```

## Arguments

- model1:

  Stage 1 model

- formula:

  Formula specifying observed covariates in stage 2 model

- data:

  data.frame

- predictfun:

  Predictions to be used in stage 2

- id1:

  Optional id-vector of stage 1

- id2:

  Optional id-vector of stage 2

- ...:

  Additional arguments to lower level functions

## See also

stack.estimate

## Examples

``` r
m <- lvm(c(y1,y2,y3)~u,c(y3,y4,y5)~v,u~~v,c(u,v)~x)
transform(m,u2~u) <- function(x) x^2
transform(m,uv~u+v) <- prod
regression(m) <- z~u2+u+v+uv+x
set.seed(1)
d <- sim(m,1000,p=c("u,u"=1))

## Stage 1
m1 <- lvm(c(y1[0:s],y2[0:s],y3[0:s])~1*u,c(y3[0:s],y4[0:s],y5[0:s])~1*v,u~b*x,u~~v)
latent(m1) <- ~u+v
e1 <- estimate(m1,d)

pp <- function(mu,var,data,...) {
    cbind(u=mu[,"u"],u2=mu[,"u"]^2+var["u","u"],
          v=mu[,"v"],uv=mu[,"u"]*mu[,"v"]+var["u","v"])
}
(e <- measurement.error(e1, z~1+x, data=d, predictfun=pp))
#>             Estimate Std.Err     2.5% 97.5%   P-value
#> (Intercept)   0.1358  0.1185 -0.09636 0.368 2.516e-01
#> x             1.1287  0.1181  0.89722 1.360 1.210e-21
#> u             0.9437  0.1217  0.70518 1.182 8.900e-15
#> u2            0.9374  0.0965  0.74830 1.127 2.615e-22
#> v             1.1385  0.1010  0.94056 1.336 1.746e-29
#> uv            1.0375  0.1079  0.82608 1.249 6.807e-22

## uu <- seq(-1,1,length.out=100)
## pp <- estimate(e,function(p,...) p["(Intercept)"]+p["u"]*uu+p["u2"]*uu^2)$coefmat
if (interactive()) {
    plot(e,intercept=TRUE,line=0)
    dev.off()
    f <- function(p) p[1]+p["u"]*u+p["u2"]*u^2
    u <- seq(-1,1,length.out=100)
    plot(e, f, data=data.frame(u), ylim=c(-.5,2.5))
}
```
