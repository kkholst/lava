# Concatenate estimate objects

When all arguments are `estimate` objects, they are merged into a single
`estimate` object. When some arguments are not `estimate` objects but
are named numeric scalars/vectors, the merged `estimate` object is
returned with an `"extra"` attribute containing the provided values.
This is useful for bundling auxiliary per-iteration information (e.g.,
convergence status) alongside an estimate for use with
[`sim.default()`](https://kkholst.github.io/lava/reference/sim.default.md).

## Usage

``` r
# S3 method for class 'estimate'
c(..., as.list = FALSE)
```

## Arguments

- ...:

  `estimate` objects and/or named numeric values

- as.list:

  if TRUE the returned object will be of class `list` and not an
  `estimate` object.

## Value

An `estimate` object. If extra named numeric values were provided, the
result carries an `"extra"` attribute (a named numeric vector).

## Details

arguments `drop.ic`, `paired`, `sep` are passed to
[merge.estimate](https://kkholst.github.io/lava/reference/merge.estimate.md)

## See also

[`sim.default()`](https://kkholst.github.io/lava/reference/sim.default.md)
[`merge.estimate()`](https://kkholst.github.io/lava/reference/merge.estimate.md)

## Examples

``` r
e <- estimate(coef = c(a = 1, b = 2), vcov = diag(2) * 0.1)
# Bundle estimate with extra information
res <- c(e, converged = 1, niter = 10)
attr(res, "extra")
#> converged     niter 
#>         1        10 
```
