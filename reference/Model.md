# Extract model

Extract or replace model object

## Usage

``` r
Model(x, ...)

Model(x, ...) <- value
```

## Arguments

- x:

  Fitted model

- ...:

  Additional arguments to be passed to the low level functions

- value:

  New model object (e.g. `lvm` or `multigroup`)

## Value

Returns a model object (e.g. `lvm` or `multigroup`)

## See also

[`Graph`](https://kkholst.github.io/lava/reference/Graph.md)

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm(y~x)
e <- estimate(m, sim(m,100))
Model(e)
#> Latent Variable Model
#>                   
#>   y ~ x   gaussian
#> 
#> Exogenous variables:                   
#>   x        gaussian
#> 
```
