# Summary method for 'sim' objects

Summary method for 'sim' objects

## Usage

``` r
# S3 method for class 'sim'
summary(
  object,
  estimate = NULL,
  se = NULL,
  confint,
  true = NULL,
  fun,
  names = NULL,
  unique.names = TRUE,
  minimal = FALSE,
  level = 0.95,
  quantiles = c(0, 0.025, 0.5, 0.975, 1),
  df = Inf,
  ...
)
```

## Arguments

- object:

  sim object

- estimate:

  (optional) columns with estimates

- se:

  (optional) columns with standard error estimates

- confint:

  (optional) list of pairs of columns with confidence limits

- true:

  (optional) vector of true parameter values

- fun:

  (optional) summary function

- names:

  (optional) names of estimates

- unique.names:

  if TRUE, unique.names will be applied to column names

- minimal:

  if TRUE, minimal summary will be returned

- level:

  confidence level (0.95)

- quantiles:

  quantiles (0,0.025,0.5,0.975,1)

- df:

  degrees of freedom in t-distribution used for constructing CIs
  (default Gaussian approximation)

- ...:

  additional levels to lower-level functions
