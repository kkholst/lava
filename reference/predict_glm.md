# Predict from a GLM with modified coefficients

Compute predictions from a fitted
[stats::glm](https://rdrr.io/r/stats/glm.html) object, optionally
substituting new parameter values. Unlike
[stats::predict.glm](https://rdrr.io/r/stats/predict.glm.html), this
function allows the user to supply an arbitrary coefficient vector `p`,
which is useful for computing predictions at counterfactual parameter
values (e.g., during optimization or simulation). The returned value
also includes a `"grad"` attribute containing the Jacobian of
predictions with respect to the coefficients.

## Usage

``` r
predict_glm(
  x,
  p = coef(x),
  data,
  offset = NULL,
  type = c("response", "link"),
  ...
)
```

## Arguments

- x:

  A fitted `glm` object.

- p:

  Numeric vector of coefficients (defaults to `coef(x)`).

- data:

  Optional data frame for computing the model matrix and response. If
  missing, the original model matrix from the fitted object is used.

- offset:

  Optional offset vector. If `NULL` (default), the offset stored in the
  fitted model (`x$offset`) is used.

- type:

  Character; `"response"` (default) returns predictions on the response
  scale via the inverse link function, `"link"` returns predictions on
  the linear predictor scale.

- ...:

  Additional arguments (currently unused).

## Value

Numeric vector of predictions with a `"grad"` attribute containing the
gradient (Jacobian) of predictions with respect to the coefficients.
When `type = "link"`, the gradient is simply the model matrix `X`.

## See also

[stats::predict.glm](https://rdrr.io/r/stats/predict.glm.html),
[score.glm](https://kkholst.github.io/lava/reference/internal.md)

## Examples

``` r
m <- glm(mpg ~ hp + wt, data = mtcars)
p0 <- coef(m)
# Predictions at fitted coefficients match predict.glm
all.equal(as.numeric(predict_glm(m)), fitted(m))
#> [1] "names for current but not for target"
# Predictions with modified coefficients
predict_glm(m, p = p0 * 1.1)
#>                         [,1]
#> Mazda RX4           25.92956
#> Mazda RX4 Wag       24.84183
#> Datsun 710          27.80340
#> Hornet 4 Drive      23.39152
#> Hornet Sportabout   20.15999
#> Valiant             22.52120
#> Duster 360          17.15895
#> Merc 240D           25.17577
#> Merc 230            24.19304
#> Merc 280            21.97741
#> Merc 280C           21.97741
#> Merc 450SE          17.29791
#> Merc 450SL          18.74821
#> Merc 450SLC         18.53493
#> Cadillac Fleetwood  11.39073
#> Lincoln Continental 10.29901
#> Chrysler Imperial   10.11174
#> Fiat 128            29.25893
#> Honda Civic         32.24362
#> Toyota Corolla      30.85083
#> Toyota Corona       27.04509
#> Dodge Challenger    20.69250
#> AMC Javelin         21.05508
#> Camaro Z28          16.00723
#> Pontiac Firebird    18.43242
#> Fiat X1-9           30.38932
#> Porsche 914-2       28.64111
#> Lotus Europa        30.54675
#> Ford Pantera L      18.20114
#> Ferrari Dino        23.01795
#> Maserati Bora       14.01342
#> Volvo 142E          25.28201
#> attr(,"grad")
#>                     (Intercept)  hp    wt
#> Mazda RX4                     1 110 2.620
#> Mazda RX4 Wag                 1 110 2.875
#> Datsun 710                    1  93 2.320
#> Hornet 4 Drive                1 110 3.215
#> Hornet Sportabout             1 175 3.440
#> Valiant                       1 105 3.460
#> Duster 360                    1 245 3.570
#> Merc 240D                     1  62 3.190
#> Merc 230                      1  95 3.150
#> Merc 280                      1 123 3.440
#> Merc 280C                     1 123 3.440
#> Merc 450SE                    1 180 4.070
#> Merc 450SL                    1 180 3.730
#> Merc 450SLC                   1 180 3.780
#> Cadillac Fleetwood            1 205 5.250
#> Lincoln Continental           1 215 5.424
#> Chrysler Imperial             1 230 5.345
#> Fiat 128                      1  66 2.200
#> Honda Civic                   1  52 1.615
#> Toyota Corolla                1  65 1.835
#> Toyota Corona                 1  97 2.465
#> Dodge Challenger              1 150 3.520
#> AMC Javelin                   1 150 3.435
#> Camaro Z28                    1 245 3.840
#> Pontiac Firebird              1 175 3.845
#> Fiat X1-9                     1  66 1.935
#> Porsche 914-2                 1  91 2.140
#> Lotus Europa                  1 113 1.513
#> Ford Pantera L                1 264 3.170
#> Ferrari Dino                  1 175 2.770
#> Maserati Bora                 1 335 3.570
#> Volvo 142E                    1 109 2.780
```
