# Add non-linear constraints to latent variable model

Add non-linear constraints to latent variable model

## Usage

``` r
# Default S3 method
constrain(x, par, args, endogenous = TRUE, ...) <- value

# S3 method for class 'multigroup'
constrain(x, par, k = 1, ...) <- value

constraints(object,data=model.frame(object),vcov=object$vcov,level=0.95,
                        p=pars.default(object),k,idx,...)
```

## Arguments

- x:

  `lvm`-object

- ...:

  Additional arguments to be passed to the low level functions

- value:

  Real function taking args as a vector argument

- par:

  Name of new parameter. Alternatively a formula with lhs specifying the
  new parameter and the rhs defining the names of the parameters or
  variable names defining the new parameter (overruling the `args`
  argument).

- args:

  Vector of variables names or parameter names that are used in defining
  `par`

- endogenous:

  TRUE if variable is endogenous (sink node)

- k:

  For multigroup models this argument specifies which group to
  add/extract the constraint

- object:

  `lvm`-object

- data:

  Data-row from which possible non-linear constraints should be
  calculated

- vcov:

  Variance matrix of parameter estimates

- level:

  Level of confidence limits

- p:

  Parameter vector

- idx:

  Index indicating which constraints to extract

## Value

A `lvm` object.

## Details

Add non-linear parameter constraints as well as non-linear associations
between covariates and latent or observed variables in the model
(non-linear regression).

As an example we will specify the follow multiple regression model:

\$\$E(Y\|X_1,X_2) = \alpha + \beta_1 X_1 + \beta_2 X_2\$\$
\$\$V(Y\|X_1,X_2) = v\$\$

which is defined (with the appropiate parameter labels) as

`m <- lvm(y ~ f(x,beta1) + f(x,beta2))`

`intercept(m) <- y ~ f(alpha)`

`covariance(m) <- y ~ f(v)`

The somewhat strained parameter constraint \$\$ v =
\frac{(beta1-beta2)^2}{alpha}\$\$

can then specified as

`constrain(m,v ~ beta1 + beta2 + alpha) <- function(x) (x[1]-x[2])^2/x[3] `

A subset of the arguments `args` can be covariates in the model,
allowing the specification of non-linear regression models. As an
example the non-linear regression model \$\$ E(Y\mid X) = \nu +
\Phi(\alpha + \beta X)\$\$ where \\\Phi\\ denotes the standard normal
cumulative distribution function, can be defined as

`m <- lvm(y ~ f(x,0)) # No linear effect of x`

Next we add three new parameters using the `parameter` assigment
function:

`parameter(m) <- ~nu+alpha+beta`

The intercept of \\Y\\ is defined as `mu`

`intercept(m) <- y ~ f(mu)`

And finally the newly added intercept parameter `mu` is defined as the
appropiate non-linear function of \\\alpha\\, \\\nu\\ and \\\beta\\:

`constrain(m, mu ~ x + alpha + nu) <- function(x) pnorm(x[1]*x[2])+x[3]`

The `constraints` function can be used to show the estimated non-linear
parameter constraints of an estimated model object (`lvmfit` or
`multigroupfit`). Calling `constrain` with no additional arguments
beyound `x` will return a list of the functions and parameter names
defining the non-linear restrictions.

The gradient function can optionally be added as an attribute `grad` to
the return value of the function defined by `value`. In this case the
analytical derivatives will be calculated via the chain rule when
evaluating the corresponding score function of the log-likelihood. If
the gradient attribute is omitted the chain rule will be applied on a
numeric approximation of the gradient.

## See also

[`regression`](https://kkholst.github.io/lava/reference/regression-set.md),
[`intercept`](https://kkholst.github.io/lava/reference/intercept.md),
[`covariance`](https://kkholst.github.io/lava/reference/covariance.md)

## Author

Klaus K. Holst

## Examples

``` r
## Non-linear parameter constraints 1
m <- lvm(y ~ f(x1,gamma)+f(x2,beta))
covariance(m) <- y ~ f(v)
d <- sim(m,100)
m1 <- m; constrain(m1,beta ~ v) <- function(x) x^2
## Define slope of x2 to be the square of the residual variance of y
## Estimate both restricted and unrestricted model
e <- estimate(m,d,control=list(method="NR"))
e1 <- estimate(m1,d)
p1 <- coef(e1)
p1 <- c(p1[1:2],p1[3]^2,p1[3])
## Likelihood of unrestricted model evaluated in MLE of restricted model
logLik(e,p1)
#> 'log Lik.' -144.784 (df=4)
## Likelihood of restricted model (MLE)
logLik(e1)
#> 'log Lik.' -144.784 (df=3)

## Non-linear regression

## Simulate data
m <- lvm(c(y1,y2)~f(x,0)+f(eta,1))
latent(m) <- ~eta
covariance(m,~y1+y2) <- "v"
intercept(m,~y1+y2) <- "mu"
covariance(m,~eta) <- "zeta"
intercept(m,~eta) <- 0
set.seed(1)
d <- sim(m,100,p=c(v=0.01,zeta=0.01))[,manifest(m)]
d <- transform(d,
               y1=y1+2*pnorm(2*x),
               y2=y2+2*pnorm(2*x))

## Specify model and estimate parameters
constrain(m, mu ~ x + alpha + nu + gamma) <- function(x) x[4]*pnorm(x[3]+x[1]*x[2])
if (FALSE)  ## Reduce Ex.Timings
e <- estimate(m,d,control=list(trace=1,constrain=TRUE))
constraints(e,data=d)
#> NULL
## Plot model-fit
plot(y1~x,d,pch=16); points(y2~x,d,pch=16,col="gray")
x0 <- seq(-4,4,length.out=100)
lines(x0,coef(e)["nu"] + coef(e)["gamma"]*pnorm(coef(e)["alpha"]*x0))

 # \dontrun{}

## Multigroup model
m1 <- lvm(y ~ f(x,beta)+f(z,beta2))
m2 <- lvm(y ~ f(x,psi) + z)
# simulate data from the two models
d1 <- sim(m1,500)
d2 <- sim(m2,500)
# Add 'non'-linear parameter constraint
constrain(m2,psi ~ beta2) <- function(x) x
## Add parameter beta2 to model 2, now beta2 exists in both models
parameter(m2) <- ~ beta2
ee <- estimate(list(m1,m2),list(d1,d2),control=list(method="NR"))
summary(ee)
#> ||score||^2= 2.641697e-17 
#> Latent variables: 
#> ____________________________________________________
#> Group 1 (n=500)
#>                     Estimate Std. Error  Z value Pr(>|z|)
#> Regressions:                                             
#>    y~x               0.93621    0.04767 19.64130   <1e-12
#>    y~z               1.03489    0.03217 32.17115   <1e-12
#> Intercepts:                                              
#>    y                -0.05578    0.04814 -1.15869   0.2466
#> Residual Variances:                                      
#>    y                 1.15827    0.07326 15.81139         
#> ____________________________________________________
#> Group 2 (n=500)
#>                        Estimate Std. Error  Z value Pr(>|z|)
#> Regressions:                                                
#>    y~z                  0.96338    0.04278 22.51997   <1e-12
#> Intercepts:                                                 
#>    y                   -0.01359    0.04657 -0.29191   0.7704
#> Additional Parameters:                                      
#>    beta2                1.03489    0.03217 32.17115   <1e-12
#> Residual Variances:                                         
#>    y                    1.07306    0.06787 15.81139         
#> 
#> ────────────────────────────────────────────────────────────────────────────────
#> Non-linear constraints:
#>     Estimate Std. Error  Z value      Pr(>|z|)      2.5%    97.5%
#> psi 1.034886 0.03216812 32.17115 4.470275e-227 0.9718373 1.097934
#> ────────────────────────────────────────────────────────────────────────────────
#> Estimator: gaussian 
#> ────────────────────────────────────────────────────────────────────────────────
#> 
#>  Number of observations = 1000 
#>  BIC = 2994.955 
#>  AIC = 2960.6 
#>  log-Likelihood of model = -1473.3 
#> 
#>  log-Likelihood of saturated model = -1473.128 
#>  Chi-squared statistic: q = 0.3442204 , df = 1 
#>   P(Q>q) = 0.5574032 
#> 
#>  RMSEA (90% CI): 0 (0;0.0698)
#>   P(RMSEA<0.05)=0.85508
#> 
#> rank(Information) = 7 (p=7)
#> condition(Information) = 5.215445
#> mean(score^2) = 3.773853e-18 
#> ────────────────────────────────────────────────────────────────────────────────

m3 <- lvm(y ~ f(x,beta)+f(z,beta2))
m4 <- lvm(y ~ f(x,beta2) + z)
e2 <- estimate(list(m3,m4),list(d1,d2),control=list(method="NR"))
e2
#> ____________________________________________________
#> Group 1 (n=500)
#>                     Estimate Std. Error  Z value Pr(>|z|)
#> Regressions:                                             
#>    y~x               0.93621    0.04767 19.64130   <1e-12
#>    y~z               1.03489    0.03217 32.17115   <1e-12
#> Intercepts:                                              
#>    y                -0.05578    0.04814 -1.15869   0.2466
#> Residual Variances:                                      
#>    y                 1.15827    0.07326 15.81139         
#> ____________________________________________________
#> Group 2 (n=500)
#>                     Estimate Std. Error  Z value Pr(>|z|)
#> Regressions:                                             
#>    y~x               1.03489    0.03217 32.17115   <1e-12
#>    y~z               0.96338    0.04278 22.51997   <1e-12
#> Intercepts:                                              
#>    y                -0.01359    0.04657 -0.29191   0.7704
#> Residual Variances:                                      
#>    y                 1.07306    0.06787 15.81139         
#> 
```
