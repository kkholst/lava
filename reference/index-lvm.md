# Extract the parameter indicies of a lvm object

Extracts the matrices with indices of model parameters from a latent
variable model (`lvm`). Returns a list with

- A: Matrix with fixed parameters and ones where parameters are free

- J: Manifest variable selection matrix

- M0: Index of free regression parameters

- M1: Index of free and *unique* regression parameters

- P: Matrix with fixed variance parameters and ones where parameters are
  free

- P0: Index of free variance parameters

- P1: Index of free and *unique* regression parameters

- npar.var: Number of covariance parameters

## Usage

``` r
# S3 method for class 'lvm'
index(x, ...)
```

## Arguments

- x:

  object on which to set the index

- ...:

  further arguments to be passed to or from other methods.

- value:

  new index

## See also

[`modelPar()`](https://kkholst.github.io/lava/reference/internal.md)
