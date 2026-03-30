# Regression model for binomial data with unkown group of immortals

Regression model for binomial data with unkown group of immortals
(zero-inflated binomial regression)

## Usage

``` r
zibreg(
  formula,
  formula.p = ~1,
  data,
  family = stats::binomial(),
  offset = NULL,
  start,
  var = "hessian",
  ...
)
```

## Arguments

- formula:

  Formula specifying

- formula.p:

  Formula for model of disease prevalence

- data:

  data frame

- family:

  Distribution family (see the help page `family`)

- offset:

  Optional offset

- start:

  Optional starting values

- var:

  Type of variance (robust, expected, hessian, outer)

- ...:

  Additional arguments to lower level functions

## Author

Klaus K. Holst

## Examples

``` r
## Simulation
n <- 2e3
x <- runif(n,0,20)
age <- runif(n,10,30)
z0 <- rnorm(n,mean=-1+0.05*age)
z <- cut(z0,breaks=c(-Inf,-1,0,1,Inf))
p0 <- lava:::expit(model.matrix(~z+age) %*% c(-.4, -.4, 0.2, 2, -0.05))
y <- (runif(n)<lava:::tigol(-1+0.25*x-0*age))*1
u <- runif(n)<p0
y[u==0] <- 0
d <- data.frame(y=y,x=x,u=u*1,z=z,age=age)
head(d)
#>   y         x u         z      age
#> 1 0 10.616487 0    (-1,0] 27.67200
#> 2 1 18.854095 1 (-Inf,-1] 12.28950
#> 3 0 14.244491 0    (-1,0] 16.48295
#> 4 0 14.489812 0    (-1,0] 17.33514
#> 5 0  9.402572 0    (-1,0] 15.07897
#> 6 0  2.405645 0  (1, Inf] 12.41499

## Estimation
e0 <- zibreg(y~x*z,~1+z+age,data=d)
e <- zibreg(y~x,~1+z+age,data=d)
compare(e,e0)
#> 
#>  - Likelihood ratio test -
#> 
#> data:  
#> chisq = 3.8978, df = 6, p-value = 0.6905
#> sample estimates:
#> log likelihood (model 1) log likelihood (model 2) 
#>                -857.7595                -855.8106 
#> 
e
#>                   Estimate       2.5%       97.5%      P-value
#> (Intercept)    -1.04506702 -1.6056511 -0.48448289 2.583311e-04
#> x               0.26955285  0.0945194  0.44458631 2.541473e-03
#> pr:(Intercept) -0.03873774 -0.6427755  0.56530001 8.999733e-01
#> pr:z(-1,0]     -0.39644128 -0.8162485  0.02336592 6.418754e-02
#> pr:z(0,1]       0.50741986  0.1012374  0.91360231 1.434653e-02
#> pr:z(1, Inf]    2.33022119  1.7820864  2.87835593 7.937660e-17
#> pr:age         -0.07884551 -0.1048220 -0.05286903 2.697676e-09
#> 
#> Prevalence probabilities:
#>                              Estimate      2.5%     97.5%
#> {(Intercept)}               0.4903168 0.3446194 0.6376780
#> {(Intercept)} + {z(-1,0]}   0.3928903 0.2637622 0.5389577
#> {(Intercept)} + {z(0,1]}    0.6150718 0.4548615 0.7536939
#> {(Intercept)} + {z(1, Inf]} 0.9081692 0.8035858 0.9598482
#> {(Intercept)} + {age}       0.4706380 0.3314733 0.6145224
PD(e0,intercept=c(1,3),slope=c(2,6))
#>     Estimate Std.Err      2.5%    97.5%
#> 50% 7.426706 5.87661 -4.091238 18.94465
#> attr(,"b")
#> [1] -1.3091313  0.1762735

B <- rbind(c(1,0,0,0,20),
           c(1,1,0,0,20),
           c(1,0,1,0,20),
           c(1,0,0,1,20))
prev <- summary(e,pr.contrast=B)$prevalence

x <- seq(0,100,length.out=100)
newdata <- expand.grid(x=x,age=20,z=levels(d$z))
fit <- predict(e,newdata=newdata)
plot(0,0,type="n",xlim=c(0,101),ylim=c(0,1),xlab="x",ylab="Probability(Event)")
count <- 0
for (i in levels(newdata$z)) {
  count <- count+1
  lines(x,fit[which(newdata$z==i)],col="darkblue",lty=count)
}
abline(h=prev[3:4,1],lty=3:4,col="gray")
abline(h=prev[3:4,2],lty=3:4,col="lightgray")
abline(h=prev[3:4,3],lty=3:4,col="lightgray")
legend("topleft",levels(d$z),col="darkblue",lty=seq_len(length(levels(d$z))))
```
