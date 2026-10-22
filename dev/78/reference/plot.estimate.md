# Plot method for 'estimate' objects

Plot method for 'estimate' objects

## Usage

``` r
# S3 method for class 'estimate'
plot(
  x,
  f,
  idx,
  intercept = FALSE,
  data,
  confint = TRUE,
  type = "l",
  xlab = "x",
  ylab = "f(x)",
  col = 1,
  add = FALSE,
  null = 0,
  ...
)
```

## Arguments

- x:

  estimate object

- f:

  function of parameter coefficients and data parsed on to 'estimate'.
  If omitted a forest-plot will be produced.

- idx:

  Index of parameters (default all)

- intercept:

  include intercept in forest-plot

- data:

  data.frame

- confint:

  Add confidence limits

- type:

  plot type ('l')

- xlab:

  x-axis label

- ylab:

  y-axis label

- col:

  color

- add:

  add plot to current device

- null:

  null value for forest-plot

- ...:

  additional arguments to lower-level functions
