# Converts strings to formula

Converts a vector of predictors and a vector of responses (characters)
i#nto a formula expression.

## Usage

``` r
toformula(y = ".", x = ".")
```

## Arguments

- y:

  vector of predictors

- x:

  vector of responses

## Value

An object of class `formula`

## See also

[`as.formula`](https://rdrr.io/r/stats/formula.html),

## Author

Klaus K. Holst

## Examples

``` r
toformula(c("age","gender"), "weight")
#> c(age, gender) ~ weight
#> <environment: 0x55809cda8768>
```
