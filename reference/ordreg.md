# Univariate cumulative link regression models

Ordinal regression models

## Usage

``` r
ordreg(
  formula,
  data = parent.frame(),
  offset,
  family = stats::binomial("probit"),
  start,
  fast = FALSE,
  ...
)
```

## Arguments

- formula:

  formula

- data:

  data.frame

- offset:

  offset

- family:

  family (default proportional odds)

- start:

  optional starting values

- fast:

  If TRUE standard errors etc. will not be calculated

- ...:

  Additional arguments to lower level functions

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm(y~x)
ordinal(m,K=3) <- ~y
d <- sim(m,100)
e <- ordreg(y~x,d)
```
