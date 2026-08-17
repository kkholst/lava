# Generic print method

Nicer print method for tabular data. Falls back to standard print method
for all other data types.

## Usage

``` r
Print(x, n = 5, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- x:

  object to print

- n:

  number of rows to show from top and bottom of tabular data

- digits:

  precision

- ...:

  additional arguments to print method
