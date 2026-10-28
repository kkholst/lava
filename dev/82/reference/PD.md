# Dose response calculation for binomial regression models

Dose response calculation for binomial regression models

## Usage

``` r
PD(
  model,
  intercept = 1,
  slope = 2,
  prob = NULL,
  x,
  level = 0.5,
  ci.level = 0.95,
  vcov,
  family,
  EB = NULL
)
```

## Arguments

- model:

  Model object or vector of parameter estimates

- intercept:

  Index of intercept parameters

- slope:

  Index of intercept parameters

- prob:

  Index of mixture parameters (only relevant for `zibreg` models)

- x:

  Optional weights
  length(x)=length(intercept)+length(slope)+length(prob)

- level:

  Probability at which level to calculate dose

- ci.level:

  Level of confidence limits

- vcov:

  Optional estimate of variance matrix of parameter estimates

- family:

  Optional distributional family argument

- EB:

  Optional ratio of treatment effect and adverse effects used to find
  optimal dose (regret-function argument)

## Author

Klaus K. Holst
