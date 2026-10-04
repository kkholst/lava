# Add an observed event time outcome to a latent variable model.

For example, if the model 'm' includes latent event time variables are
called 'T1' and 'T2' and 'C' is the end of follow-up (right censored),
then one can specify

`eventTime(object=m,formula=ObsTime~min(T1=a,T2=b,C=0,"ObsEvent"))`

when data are simulated from the model one gets 2 new columns:

- "ObsTime": the smallest of T1, T2 and C

- "ObsEvent": 'a' if T1 is smallest, 'b' if T2 is smallest and '0' if C
  is smallest

Note that "ObsEvent" and "ObsTime" are names specified by the user.

## Usage

``` r
eventTime(object, formula, eventName = "status", ...)
```

## Arguments

- object:

  Model object

- formula:

  Formula (see details)

- eventName:

  Event names

- ...:

  Additional arguments to lower levels functions

## Author

Thomas A. Gerds, Klaus K. Holst

## Examples

``` r

# Right censored survival data without covariates
m0 <- lvm()
distribution(m0,"eventtime") <- dist_cox_weibull(scale=1/100,shape=2)
distribution(m0,"censtime") <- dist_cox_exponential(rate=1/10)
m0 <- eventTime(m0,time~min(eventtime=1,censtime=0),"status")
sim(m0,10)
#>    eventtime    censtime        time status
#> 1   8.846923 18.79067364  8.84692274      1
#> 2  14.009893  3.39319492  3.39319492      0
#> 3  14.029624  4.81329027  4.81329027      0
#> 4  17.929567 22.09611099 17.92956695      1
#> 5  10.537679 16.55291245 10.53767868      1
#> 6   8.534200  0.06189218  0.06189218      0
#> 7   8.161939 11.85611149  8.16193939      1
#> 8   7.802480  0.31966069  0.31966069      0
#> 9  13.206909  6.78797323  6.78797323      0
#> 10 14.468739  2.16618630  2.16618630      0

# Alternative specification of the right censored survival outcome
## eventTime(m,"Status") <- ~min(eventtime=1,censtime=0)

# Cox regression:
# lava implements two different parametrizations of the same
# Weibull regression model. The first specifies
# the effects of covariates as proportional hazard ratios
# and works as follows:
m <- lvm()
distribution(m,"eventtime") <- dist_cox_weibull(scale=1/100,shape=2)
distribution(m,"censtime") <- dist_cox_weibull(scale=1/100,shape=2)
m <- eventTime(m,time~min(eventtime=1,censtime=0),"status")
distribution(m,"sex") <- dist_bernoulli(p=0.4)
distribution(m,"sbp") <- dist_gaussian(mean=120,sd=20)
regression(m,from="sex",to="eventtime") <- 0.4
regression(m,from="sbp",to="eventtime") <- -0.01
sim(m,6)
#>   eventtime   censtime       time status sex      sbp
#> 1  4.503515 10.3321280  4.5035153      1   0 106.6415
#> 2  7.788162  8.9109949  7.7881625      1   0 103.8519
#> 3  9.875467  2.0632497  2.0632497      0   1 134.3462
#> 4 12.386353 12.8740043 12.3863532      1   1 107.2798
#> 5 16.413136  0.4363755  0.4363755      0   1 113.6455
#> 6 18.454376 13.7010274 13.7010274      0   0 122.3273
# The parameters can be recovered using a Cox regression
# routine or a Weibull regression model. E.g.,
if (FALSE) { # \dontrun{
    set.seed(18)
    d <- sim(m,1000)
    library(survival)
    coxph(Surv(time,status)~sex+sbp,data=d)

    sr <- survreg(Surv(time,status)~sex+sbp,data=d)
    library(SurvRegCensCov)
    ConvertWeibull(sr)

} # }

# The second parametrization is an accelerated failure time
# regression model and uses the function dist_weibull instead
# of dist_cox_weibull to specify the event time distributions.
# Here is an example:

ma <- lvm()
distribution(ma,"eventtime") <- dist_weibull(scale=3,shape=1/0.7)
distribution(ma,"censtime") <- dist_weibull(scale=2,shape=1/0.7)
ma <- eventTime(ma,time~min(eventtime=1,censtime=0),"status")
distribution(ma,"sex") <- dist_bernoulli(p=0.4)
distribution(ma,"sbp") <- dist_gaussian(mean=120,sd=20)
regression(ma,from="sex",to="eventtime") <- 0.7
regression(ma,from="sbp",to="eventtime") <- -0.008
set.seed(17)
sim(ma,6)
#>   eventtime  censtime      time status sex       sbp
#> 1 0.5531481 1.1285503 0.5531481      1   1  99.69983
#> 2 4.2973225 1.4665922 1.4665922      0   1 118.40727
#> 3 1.5884110 0.4704796 0.4704796      0   1 115.34026
#> 4 1.7404946 1.2284359 1.2284359      0   1 103.65464
#> 5 0.2765550 0.8633771 0.2765550      1   1 135.44182
#> 6 1.5803203 0.6912997 0.6912997      0   0 116.68776
# The regression coefficients of the AFT model
# can be tranformed into log(hazard ratios):
#  coef.coxWeibull = - coef.weibull / shape.weibull
if (FALSE) { # \dontrun{
    set.seed(17)
    da <- sim(ma,1000)
    library(survival)
    fa <- coxph(Surv(time,status)~sex+sbp,data=da)
    coef(fa)
    c(0.7,-0.008)/0.7
} # }


# The following are equivalent parametrizations
# which produce exactly the same random numbers:

model.aft <- lvm()
distribution(model.aft,"eventtime") <-
  dist_weibull(intercept=-log(1/100)/2,sigma=1/2)
distribution(model.aft,"censtime") <-
  dist_weibull(intercept=-log(1/100)/2,sigma=1/2)
sim(model.aft,6,seed=17)
#>   eventtime  censtime
#> 1 12.552208 13.652847
#> 2 12.946401  1.792538
#> 3  4.984980  8.710482
#> 4 12.806975  5.025406
#> 5  9.133336  9.469785
#> 6 24.669793  7.863944

model.aft <- lvm()
distribution(model.aft,"eventtime") <- dist_weibull(scale=100^(1/2), shape=2)
distribution(model.aft,"censtime") <- dist_weibull(scale=100^(1/2), shape=2)
sim(model.aft,6,seed=17)
#>   eventtime  censtime
#> 1 12.552208 13.652847
#> 2 12.946401  1.792538
#> 3  4.984980  8.710482
#> 4 12.806975  5.025406
#> 5  9.133336  9.469785
#> 6 24.669793  7.863944

model.cox <- lvm()
distribution(model.cox,"eventtime") <- dist_cox_weibull(scale=1/100,shape=2)
distribution(model.cox,"censtime") <- dist_cox_weibull(scale=1/100,shape=2)
sim(model.cox,6,seed=17)
#>   eventtime  censtime
#> 1 12.552208 13.652847
#> 2 12.946401  1.792538
#> 3  4.984980  8.710482
#> 4 12.806975  5.025406
#> 5  9.133336  9.469785
#> 6 24.669793  7.863944

# The minimum of multiple latent times one of them still
# being a censoring time, yield
# right censored competing risks data

mc <- lvm()
distribution(mc,~X2) <- dist_bernoulli()
regression(mc) <- T1~f(X1,-.5)+f(X2,0.3)
regression(mc) <- T2~f(X2,0.6)
distribution(mc,~T1) <- dist_cox_weibull(scale=1/100)
distribution(mc,~T2) <- dist_cox_weibull(scale=1/100)
distribution(mc,~C) <- dist_cox_weibull(scale=1/100)
mc <- eventTime(mc,time~min(T1=1,T2=2,C=0),"event")
sim(mc,6)
#>   X2        T1          X1        T2         C      time event
#> 1  0  7.023211 -0.05517906 14.575911 11.814841  7.023211     1
#> 2  1  6.179319  0.83847112  5.138275  7.716248  5.138275     2
#> 3  1  5.890305  0.15937013  9.886258 12.038218  5.890305     1
#> 4  0 15.128422  0.62595440 13.871923 11.629342 11.629342     0
#> 5  1 12.075247  0.63358473  9.551212  2.196766  2.196766     0
#> 6  0 19.313957  0.68102765  7.433206  3.930187  3.930187     0

```
