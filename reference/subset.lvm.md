# Extract subset of latent variable model

Extract measurement models or user-specified subset of model

## Usage

``` r
# S3 method for class 'lvm'
subset(x, vars, ...)
```

## Arguments

- x:

  `lvm`-object.

- vars:

  Character vector or formula specifying variables to include in subset.

- ...:

  Additional arguments to be passed to the low level functions

## Value

A `lvm`-object.

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm(c(y1,y2)~x1+x2)
subset(m,~y1+x1)
#> Latent Variable Model
#>                     
#>   y1 ~ x1   gaussian
#> 
#> Exogenous variables:                    
#>   x1        gaussian
#> 
```
