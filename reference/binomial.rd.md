# Define constant risk difference or relative risk association for binary exposure

Set up model as defined in Richardson, Robins and Wang (2017).

## Usage

``` r
binomial.rd(
  x,
  response,
  exposure,
  target.model,
  nuisance.model,
  exposure.model = binomial.lvm(),
  ...
)
```

## Arguments

- x:

  model

- response:

  response variable (character or formula)

- exposure:

  exposure variable (character or formula)

- target.model:

  variable defining the linear predictor for the target model

- nuisance.model:

  variable defining the linear predictor for the nuisance model

- exposure.model:

  model for exposure (default binomial logit link)

- ...:

  additional arguments to lower level functions

## Examples

``` r
## ---------------------------------------------------------------
## binomial.rd: constant risk-difference model
##   P(Y=1|Z=1) - P(Y=1|Z=0) = tanh(lp)
## ---------------------------------------------------------------
m <- lvm()
regression(m) <- z ~ x
regression(m) <- lp ~ x
regression(m) <- op ~ x
intercept(m, ~lp) <- 0.4   ## constant linear predictor for RD
intercept(m, ~op) <- 0     ## odds product = exp(0) = 1
distribution(m, ~lp) <- normal.lvm(sd = 0)
distribution(m, ~op) <- normal.lvm(sd = 0)
m <- binomial.rd(m, response = "y", exposure = "z",
                 target.model = "lp", nuisance.model = "op")
set.seed(1)
d <- sim(m, n = 2000)
## Empirical risk difference should be close to tanh(0.4)
mean(d$y[d$z == 1]) - mean(d$y[d$z == 0])
#> [1] 0.2653421
tanh(0.4)
#> [1] 0.379949

## Formula interface: response ~ exposure | target | nuisance
m2 <- lvm()
regression(m2) <- z ~ x
regression(m2) <- lp ~ x
regression(m2) <- op ~ x
m2 <- binomial.rd(m2, y ~ z | lp | op)

## ---------------------------------------------------------------
## binomial.rr: constant relative-risk model
##   log(P(Y=1|Z=1) / P(Y=1|Z=0)) = lp
## ---------------------------------------------------------------
m <- lvm()
regression(m) <- z ~ x
regression(m) <- lp ~ x
regression(m) <- op ~ x
intercept(m, ~lp) <- log(1.5)   ## constant log relative-risk
intercept(m, ~op) <- 0          ## odds product = 1
distribution(m, ~lp) <- normal.lvm(sd = 0)
distribution(m, ~op) <- normal.lvm(sd = 0)
m <- binomial.rr(m, response = "y", exposure = "z",
                 target.model = "lp", nuisance.model = "op")
set.seed(1)
d <- sim(m, n = 2000)
## Empirical log-RR should be close to log(1.5)
log(mean(d$y[d$z == 1]) / mean(d$y[d$z == 0]))
#> [1] 0.4543279
log(1.5)
#> [1] 0.4054651
```
