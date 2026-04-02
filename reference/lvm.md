# Initialize new latent variable model

Function that constructs a new latent variable model object

## Usage

``` r
lvm(x = NULL, ..., latent = NULL, messages = lava.options()$messages)
```

## Arguments

- x:

  Vector of variable names. Optional but gives control of the sequence
  of appearance of the variables. The argument can be given as a
  character vector or formula, e.g. `~y1+y2` is equivalent to
  `c("y1","y2")`. Alternatively the argument can be a formula specifying
  a linear model.

- ...:

  Additional arguments to be passed to the low level functions

- latent:

  (optional) Latent variables

- messages:

  Controls what messages are printed (0: none)

## Value

Returns an object of class `lvm`.

## See also

[`regression`](https://kkholst.github.io/lava/reference/regression-set.md),
[`covariance`](https://kkholst.github.io/lava/reference/covariance.md),
[`intercept`](https://kkholst.github.io/lava/reference/intercept.md),
...

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm() # Empty model
m1 <- lvm(y~x) # Simple linear regression
m2 <- lvm(~y1+y2) # Model with two independent variables (argument)
m3 <- lvm(list(c(y1,y2,y3)~u,u~x+z)) # SEM with three items
```
