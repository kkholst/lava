# Predict function for latent variable models

Predictions of conditinoal mean and variance and calculation of jacobian
with respect to parameter vector.

## Usage

``` r
predictlvm(object, formula, p = coef(object), data = model.frame(object), ...)
```

## Arguments

- object:

  Model object

- formula:

  Formula specifying which variables to predict and which to condition
  on

- p:

  Parameter vector

- data:

  Data.frame

- ...:

  Additional arguments to lower level functions

## See also

predict.lvm

## Examples

``` r
m <- lvm(c(x1,x2,x3)~u1,u1~z,
         c(y1,y2,y3)~u2,u2~u1+z)
latent(m) <- ~u1+u2
d <- simulate(m,10,"u2,u2"=2,"u1,u1"=0.5,seed=123)
e <- estimate(m,d)

## Conditional mean given covariates
predictlvm(e,c(x1,x2)~1)$mean
#>                x1           x2
#>  [1,] -0.17634038  0.001097242
#>  [2,]  0.22409175  0.370702775
#>  [3,] -0.64578819 -0.432210919
#>  [4,]  2.17930239  2.175394823
#>  [5,]  1.38879089  1.445739518
#>  [6,] -0.52874258 -0.324175875
#>  [7,]  0.06371187  0.222669470
#>  [8,]  0.01125438  0.174250335
#>  [9,]  1.03672161  1.120773697
#> [10,]  0.32654483  0.465268679
## Conditional variance of u1,y1 given x1,x2
predictlvm(e,c(u1,y1)~x1+x2)$var
#>            u1         y1
#> u1 0.17501758 0.09920436
#> y1 0.09920436 0.89761525
```
