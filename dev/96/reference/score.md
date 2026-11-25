# Extract score function

Extract the score function (gradient of the log-likelihood) from a model
object.

## Usage

``` r
score(x, ...)
```

## Arguments

- x:

  Model object

- ...:

  Additional arguments to lower level functions

## Value

Matrix of score/gradient values. Rows correspond to observations (if
`indiv=TRUE`) or a single row, columns to parameters.
