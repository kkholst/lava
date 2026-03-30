# Weighted K-means

Weighted K-means via Lloyd's algorithm

## Usage

``` r
wkm(
  x,
  mu,
  data,
  weights = rep(1, NROW(x)),
  iter.max = 20,
  n.start = 5,
  init = "kmpp",
  ...
)
```

## Arguments

- x:

  Data (or formula)

- mu:

  Initial centers (or number centers chosen randomly among x)

- data:

  optional data frmae

- weights:

  Optional weights

- iter.max:

  Max number of iterations

- n.start:

  Number of restarts

- init:

  method to create initial centres (default kmeans++)

- ...:

  Additional arguments to lower level functions

## Author

Klaus K. Holst
