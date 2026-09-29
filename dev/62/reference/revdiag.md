# Create/extract 'reverse'-diagonal matrix or off-diagonal elements

Create/extract 'reverse'-diagonal matrix or off-diagonal elements

## Usage

``` r
revdiag(x,...)
offdiag(x,type=0,...)

revdiag(x, ...) <- value
offdiag(x, type = 0, ...) <- value
```

## Arguments

- x:

  vector

- ...:

  additional arguments to lower level functions

- value:

  For the assignment function the values to put in the diagonal

- type:

  0: upper and lower triangular, 1: upper triangular, 2: lower
  triangular, 3: upper triangular + diagonal, 4: lower triangular +
  diagonal

## Author

Klaus K. Holst
