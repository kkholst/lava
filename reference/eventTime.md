# Add an observed event time outcome to a latent variable model.

For example, if the model 'm' includes latent event time variables are
called 'T1' and 'T2' and 'C' is the end of follow-up (right censored),
then one can specify

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

## Details

`eventTime(object=m,formula=ObsTime~min(T1=a,T2=b,C=0,"ObsEvent"))`

when data are simulated from the model one gets 2 new columns:

- "ObsTime": the smallest of T1, T2 and C

- "ObsEvent": 'a' if T1 is smallest, 'b' if T2 is smallest and '0' if C
  is smallest

Note that "ObsEvent" and "ObsTime" are names specified by the user.

## Author

Thomas A. Gerds, Klaus K. Holst

## Examples

``` r
# Right censored survival data without covariates
m0 <- lvm()
distribution(m0,"eventtime") <- coxWeibull.lvm(scale=1/100,shape=2)
distribution(m0,"censtime") <- coxExponential.lvm(rate=1/10)
m0 <- eventTime(m0,time~min(eventtime=1,censtime=0),"status")
sim(m0,10)
#>    eventtime   censtime       time status
#> 1   7.775310 13.7796884  7.7753101      1
#> 2  12.285109  7.0934129  7.0934129      0
#> 3  15.396895 15.6199780 15.3968954      1
#> 4  15.795398  4.7664010  4.7664010      0
#> 5   5.351330  2.6996180  2.6996180      0
#> 6   2.158865  8.0832631  2.1588648      1
#> 7   6.828698 13.3793979  6.8286979      1
#> 8  13.521635  1.9044909  1.9044909      0
#> 9  10.332128  9.7990463  9.7990463      0
#> 10  8.910995  0.7311118  0.7311118      0

# Alternative specification of the right censored survival outcome
## eventTime(m,"Status") <- ~min(eventtime=1,censtime=0)

# Cox regression:
# lava implements two different parametrizations of the same
# Weibull regression model. The first specifies
# the effects of covariates as proportional hazard ratios
# and works as follows:
m <- lvm()
distribution(m,"eventtime") <- coxWeibull.lvm(scale=1/100,shape=2)
distribution(m,"censtime") <- coxWeibull.lvm(scale=1/100,shape=2)
m <- eventTime(m,time~min(eventtime=1,censtime=0),"status")
distribution(m,"sex") <- binomial.lvm(p=0.4)
distribution(m,"sbp") <- normal.lvm(mean=120,sd=20)
regression(m,from="sex",to="eventtime") <- 0.4
regression(m,from="sbp",to="eventtime") <- -0.01
sim(m,6)
#>   eventtime  censtime      time status sex      sbp
#> 1 36.921926  7.448510  7.448510      0   0 154.6311
#> 2 18.130309  6.757759  6.757759      0   1 177.8781
#> 3 35.747975  2.942973  2.942973      0   1 149.9032
#> 4  8.547178 14.147578  8.547178      1   1 129.5842
#> 5 11.795046  1.662079  1.662079      0   0 108.0632
#> 6 11.760196 20.535047 11.760196      1   0 117.8530
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
# regression model and uses the function weibull.lvm instead
# of coxWeibull.lvm to specify the event time distributions.
# Here is an example:

ma <- lvm()
distribution(ma,"eventtime") <- weibull.lvm(scale=3,shape=1/0.7)
distribution(ma,"censtime") <- weibull.lvm(scale=2,shape=1/0.7)
ma <- eventTime(ma,time~min(eventtime=1,censtime=0),"status")
distribution(ma,"sex") <- binomial.lvm(p=0.4)
distribution(ma,"sbp") <- normal.lvm(mean=120,sd=20)
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
distribution(model.aft,"eventtime") <- weibull.lvm(intercept=-log(1/100)/2,sigma=1/2)
distribution(model.aft,"censtime") <- weibull.lvm(intercept=-log(1/100)/2,sigma=1/2)
sim(model.aft,6,seed=17)
#>   eventtime  censtime
#> 1 12.552208 13.652847
#> 2 12.946401  1.792538
#> 3  4.984980  8.710482
#> 4 12.806975  5.025406
#> 5  9.133336  9.469785
#> 6 24.669793  7.863944

model.aft <- lvm()
distribution(model.aft,"eventtime") <- weibull.lvm(scale=100^(1/2), shape=2)
distribution(model.aft,"censtime") <- weibull.lvm(scale=100^(1/2), shape=2)
sim(model.aft,6,seed=17)
#>   eventtime  censtime
#> 1 12.552208 13.652847
#> 2 12.946401  1.792538
#> 3  4.984980  8.710482
#> 4 12.806975  5.025406
#> 5  9.133336  9.469785
#> 6 24.669793  7.863944

model.cox <- lvm()
distribution(model.cox,"eventtime") <- coxWeibull.lvm(scale=1/100,shape=2)
distribution(model.cox,"censtime") <- coxWeibull.lvm(scale=1/100,shape=2)
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
distribution(mc,~X2) <- binomial.lvm()
regression(mc) <- T1~f(X1,-.5)+f(X2,0.3)
regression(mc) <- T2~f(X2,0.6)
distribution(mc,~T1) <- coxWeibull.lvm(scale=1/100)
distribution(mc,~T2) <- coxWeibull.lvm(scale=1/100)
distribution(mc,~C) <- coxWeibull.lvm(scale=1/100)
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
