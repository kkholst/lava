# Concatenate summary.estimate objects

When all arguments are `summary.estimate` objects, they are merged into
a single `summary.estimate` object. When some arguments are not
`summary.estimate` objects but are named numeric scalars/vectors, an
object of class `summary.estimate` is returned with an additional
`extra` attribute that contains the provided numeric scalars/vectors.
This is useful for bundling auxiliary per-iteration information (e.g.,
convergence status) alongside an estimate for use with
[`sim.default()`](https://kkholst.github.io/lava/reference/sim.default.md).

## Usage

``` r
# S3 method for class 'summary.estimate'
c(...)
```

## Arguments

- ...:

  `summary.estimate` objects and/or named numeric values

## Examples

``` r
e1 <- estimate(coef = 1, IC = scale(rnorm(10)), id = 1:10, labels = "a1")
e2 <- estimate(coef = 2, IC = scale(rnorm(10)), id = 1:10, labels = "a2")

# concatenating two summary.estimate objects
c(e1, e2)
#>    Estimate Std.Err  2.5% 97.5%   P-value
#> a1        1     0.3 0.412 1.588 8.581e-04
#> a2        2     0.3 1.412 2.588 2.617e-11

# concatenating one summary.estimate object with one numerical variable
ss <- c(summary(e1), niter = 2)
print(ss)
#> Concatenated summary.estimate objects: 
#> ────────────────────────────────────────────────────────────
#>    Estimate Std.Err  2.5% 97.5%   P-value
#> a1        1     0.3 0.412 1.588 0.0008581
#> ────────────────────────────────────────────────────────────
#> niter 
#>     2 
attributes(ss)$extra
#> niter 
#>     2 
```
