# Missing data example

Simulated data generated from model \$\$E(Y_i\mid X) = X, \quad
cov(Y_1,Y_2\mid X)=0.5\$\$

## Format

list of data.frames

## Source

Simulated

## Details

The list contains four data sets 1) Complete data 2) MCAR 3) MAR 4) MNAR
(missing mechanism depends on variable V correlated with Y1,Y2)

## Examples

``` r
data(missingdata)
e0 <- estimate(lvm(c(y1,y2)~b*x,y1~~y2),missingdata[[1]]) ## No missing
e1 <- estimate(lvm(c(y1,y2)~b*x,y1~~y2),missingdata[[2]]) ## CC (MCAR)
e2 <- estimate(lvm(c(y1,y2)~b*x,y1~~y2),missingdata[[2]],missing=TRUE) ## MCAR
e3 <- estimate(lvm(c(y1,y2)~b*x,y1~~y2),missingdata[[3]]) ## CC (MAR)
e4 <- estimate(lvm(c(y1,y2)~b*x,y1~~y2),missingdata[[3]],missing=TRUE) ## MAR
```
