# Spaghetti plot

Spaghetti plot for longitudinal data

## Usage

``` r
spaghetti(
  formula,
  data = NULL,
  id = "id",
  group = NULL,
  type = "o",
  lty = 1,
  pch = NA,
  col = 1:10,
  alpha = 0.3,
  lwd = 1,
  level = 0.95,
  trend.formula = formula,
  tau = NULL,
  trend.lty = 1,
  trend.join = TRUE,
  trend.delta = 0.2,
  trend = !is.null(tau),
  trend.col = col,
  trend.alpha = 0.2,
  trend.lwd = 3,
  trend.jitter = 0,
  legend = NULL,
  by = NULL,
  xlab = "Time",
  ylab = "",
  add = FALSE,
  ...
)
```

## Arguments

- formula:

  Formula (response ~ time)

- data:

  data.frame

- id:

  Id variable

- group:

  group variable

- type:

  Type (line 'l', stair 's', ...)

- lty:

  Line type

- pch:

  Colour

- col:

  Colour

- alpha:

  transparency (0-1)

- lwd:

  Line width

- level:

  Confidence level

- trend.formula:

  Formula for trendline

- tau:

  Quantile to estimate (trend)

- trend.lty:

  Trend line type

- trend.join:

  Trend polygon

- trend.delta:

  Length of limit bars

- trend:

  Add trend line

- trend.col:

  Colour of trend line

- trend.alpha:

  Transparency

- trend.lwd:

  Trend line width

- trend.jitter:

  Jitter amount

- legend:

  Legend

- by:

  make separate plot for each level in 'by' (formula, name of column, or
  vector)

- xlab:

  Label of X-axis

- ylab:

  Label of Y-axis

- add:

  Add to existing device

- ...:

  Additional arguments to lower level arguments

## Author

Klaus K. Holst

## Examples

``` r
if (interactive() & requireNamespace("mets")) {
K <- 5
y <- "y"%++%seq(K)
m <- lvm()
regression(m,y=y,x=~u) <- 1
regression(m,y=y,x=~s) <- seq(K)-1
regression(m,y=y,x=~x) <- "b"
N <- 50
d <- sim(m,N); d$z <- rbinom(N,1,0.5)
dd <- mets::fast.reshape(d); dd$num <- dd$num+3
spaghetti(y~num,dd,id="id",lty=1,col=Col(1,.4),
          trend.formula=~factor(num),trend=TRUE,trend.col="darkblue")
dd$num <- dd$num+rnorm(nrow(dd),sd=0.5) ## Unbalance
spaghetti(y~num,dd,id="id",lty=1,col=Col(1,.4),
          trend=TRUE,trend.col="darkblue")
spaghetti(y~num,dd,id="id",lty=1,col=Col(1,.4),
           trend.formula=~num+I(num^2),trend=TRUE,trend.col="darkblue")
}
```
