# Define variables as ordinal

Define variables as ordinal in latent variable model object

## Usage

``` r
ordinal(x, ...) <- value
```

## Arguments

- x:

  Object

- ...:

  additional arguments to lower level functions

- value:

  variable (formula or character vector)

## Examples

``` r
if (requireNamespace("mets")) {
m <- lvm(y + z ~ x + 1*u[0], latent=~u)
ordinal(m, K=3) <- ~y+z
d <- sim(m, 100, seed=1)
e <- estimate(m, d)
}
```
