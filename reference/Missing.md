# Missing value generator

Missing value generator

## Usage

``` r
Missing(object, formula, Rformula, missing.name, suffix = "0", ...)
```

## Arguments

- object:

  `lvm`-object.

- formula:

  The right hand side specifies the name of a latent variable which is
  not always observed. The left hand side specifies the name of a new
  variable which is equal to the latent variable but has missing values.
  If given as a string then this is used as the name of the latent
  (full-data) name, and the observed data name is 'missing.data'

- Rformula:

  Missing data mechanism with left hand side specifying the name of the
  observed data indicator (may also just be given as a character instead
  of a formula)

- missing.name:

  Name of observed data variable (only used if 'formula' was given as a
  character specifying the name of the full-data variable)

- suffix:

  If missing.name is missing, then the name of the oberved data variable
  will be the name of the full-data variable + the suffix

- ...:

  Passed to binomial.lvm.

## Value

lvm object

## Details

This function adds a binary variable to a given `lvm` model and also a
variable which is equal to the original variable where the binary
variable is equal to zero

## Author

Thomas A. Gerds \<tag@biostat.ku.dk\>

## Examples

``` r
library(lava)
set.seed(17)
m <- lvm(y0~x01+x02+x03)
m <- Missing(m,formula=x1~x01,Rformula=R1~0.3*x02+-0.7*x01,p=0.4)
sim(m,10)
#>            y0         x01        x02         x03 R1        x1
#> 1  -0.3307614  1.18078924  0.6810276 -1.17756957  0        NA
#> 2  -1.0786445  0.64319207 -0.6820334 -0.96016651  0        NA
#> 3  -0.5398741  1.29532187 -0.7232567 -0.87895224  0        NA
#> 4  -2.5119604  0.18791807  1.6735260 -3.55613648  1 0.1879181
#> 5   0.3507905  1.59120510 -0.5957556 -1.41674984  0        NA
#> 6   0.4902836 -0.05517906  1.1598438 -0.44876927  0        NA
#> 7   1.1528003  0.83847112  0.1174224 -0.77596771  1 0.8384711
#> 8   1.3032974  0.15937013  0.2592214 -0.83182805  0        NA
#> 9   1.3153836  0.62595440  0.3823621  0.05183012  1 0.6259544
#> 10 -0.3278672  0.63358473 -0.7114817 -0.61655131  0        NA


m <- lvm(y~1)
m <- Missing(m,"y","r")
## same as
## m <- Missing(m,y~1,r~1)
sim(m,10)
#>              y r          y0
#> 1   0.07419352 1  0.07419352
#> 2   1.75169617 0          NA
#> 3  -0.23148744 0          NA
#> 4   0.54345248 0          NA
#> 5  -0.98900140 0          NA
#> 6   0.31553146 1  0.31553146
#> 7   2.44232746 1  2.44232746
#> 8   0.54969286 1  0.54969286
#> 9  -0.02924337 1 -0.02924337
#> 10 -0.83078338 0          NA

## same as
m <- lvm(y~1)
Missing(m,"y") <- r~x
sim(m,10)
#>              y r          y0          x
#> 1   0.03054575 1  0.03054575  0.5602348
#> 2  -0.78551741 0          NA -1.7924178
#> 3   0.32544056 0          NA -1.5654169
#> 4  -0.88084355 0          NA -3.3203189
#> 5   0.20932594 0          NA  0.1547216
#> 6   0.15103295 1  0.15103295 -0.3646267
#> 7  -0.34347879 0          NA -2.4336839
#> 8   0.90587760 0          NA  0.3364643
#> 9   0.91895485 0          NA -0.6404528
#> 10 -0.55598749 1 -0.55598749  1.8211204

m <- lvm(y~1)
m <- Missing(m,"y","r",suffix=".")
## same as
## m <- Missing(m,"y","r",missing.name="y.")
## same as
## m <- Missing(m,y.~y,"r")
sim(m,10)
#>              y r          y.
#> 1   0.46345064 1  0.46345064
#> 2   0.88066627 1  0.88066627
#> 3   0.13942195 0          NA
#> 4  -0.37505483 0          NA
#> 5   0.04253903 0          NA
#> 6   0.63388029 1  0.63388029
#> 7  -1.70748846 1 -1.70748846
#> 8   0.01968655 1  0.01968655
#> 9  -0.29879637 0          NA
#> 10 -0.44176283 1 -0.44176283
```
