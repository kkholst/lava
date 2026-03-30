# Calculate partial correlations

Calculate partial correlation coefficients and confidence limits via
Fishers z-transform

## Usage

``` r
partialcor(formula, data, level = 0.95, ...)
```

## Arguments

- formula:

  formula speciying the covariates and optionally the outcomes to
  calculate partial correlation for

- data:

  data.frame

- level:

  Level of confidence limits

- ...:

  Additional arguments to lower level functions

## Value

A coefficient matrix

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm(c(y1,y2,y3)~x1+x2)
covariance(m) <- c(y1,y2,y3)~y1+y2+y3
d <- sim(m,500)
partialcor(~x1+x2,d)
#>             cor        z         pval   lowerCI   upperCI
#> y1~y2 0.5362128 13.30955 2.037164e-40 0.4704451 0.5960563
#> y1~y3 0.5179408 12.74715 3.233467e-37 0.4505088 0.5794966
#> y2~y3 0.5307336 13.13932 1.959856e-39 0.4644598 0.5910959
```
