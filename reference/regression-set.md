# Add regression association to latent variable model

Define regression association between variables in a `lvm`-object and
define linear constraints between model equations.

## Usage

``` r
# S3 method for class 'lvm'
regression(object = lvm(), to, from, fn = NA,
messages = lava.options()$messages, additive=TRUE, y, x, value, ...)
# S3 method for class 'lvm'
regression(object, to = NULL, quick = FALSE, ...) <- value
```

## Arguments

- object:

  `lvm`-object.

- ...:

  Additional arguments to be passed to the low level functions

- value:

  A formula specifying the linear constraints or if `to=NULL` a `list`
  of parameter values.

- to:

  Character vector of outcome(s) or formula object.

- from:

  Character vector of predictor(s).

- fn:

  Real function defining the functional form of predictors (for
  simulation only).

- messages:

  Controls which messages are turned on/off (0: all off)

- additive:

  If FALSE and predictor is categorical a non-additive effect is assumed

- y:

  Alias for 'to'

- x:

  Alias for 'from'

- quick:

  Faster implementation without parameter constraints

## Value

A `lvm`-object

## Details

The `regression` function is used to specify linear associations between
variables of a latent variable model, and offers formula syntax
resembling the model specification of e.g. `lm`.

For instance, to add the following linear regression model, to the
`lvm`-object, `m`: \$\$ E(Y\|X_1,X_2) = \beta_1 X_1 + \beta_2 X_2\$\$ We
can write

`regression(m) <- y ~ x1 + x2`

Multivariate models can be specified by successive calls with
`regression`, but multivariate formulas are also supported, e.g.

`regression(m) <- c(y1,y2) ~ x1 + x2`

defines \$\$ E(Y_i\|X_1,X_2) = \beta\_{1i} X_1 + \beta\_{2i} X_2 \$\$

The special function, `f`, can be used in the model specification to
specify linear constraints. E.g. to fix \\\beta_1=\beta_2\\ , we could
write

`regression(m) <- y ~ f(x1,beta) + f(x2,beta)`

The second argument of `f` can also be a number (e.g. defining an
offset) or be set to `NA` in order to clear any previously defined
linear constraints.

Alternatively, a more straight forward notation can be used:

`regression(m) <- y ~ beta*x1 + beta*x2`

All the parameter values of the linear constraints can be given as the
right handside expression of the assigment function `regression<-` (or
`regfix<-`) if the first (and possibly second) argument is defined as
well. E.g:

`regression(m,y1~x1+x2) <- list("a1","b1")`

defines \\E(Y_1\|X_1,X_2) = a1 X_1 + b1 X_2\\. The rhs argument can be a
mixture of character and numeric values (and NA's to remove
constraints).

The function `regression` (called without additional arguments) can be
used to inspect the linear constraints of a `lvm`-object.

## Note

Variables will be added to the model if not already present.

## See also

`intercept<-`, `covariance<-`, `constrain<-`, `parameter<-`, `latent<-`,
`cancel<-`, `kill<-`

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm() ## Initialize empty lvm-object
### E(y1|z,v) = beta1*z + beta2*v
regression(m) <- y1 ~ z + v
### E(y2|x,z,v) = beta*x + beta*z + 2*v + beta3*u
regression(m) <- y2 ~ f(x,beta) + f(z,beta)  + f(v,2) + u
### Clear restriction on association between y and
### fix slope coefficient of u to beta
regression(m, y2 ~ v+u) <- list(NA,"beta")

regression(m) ## Examine current linear parameter constraints
#> Regression parameters:
#>       y1 z    v y2 x    u   
#>    y1    *    *             
#>    y2    beta *    beta beta

## ## A multivariate model, E(yi|x1,x2) = beta[1i]*x1 + beta[2i]*x2:
m2 <- lvm(c(y1,y2) ~ x1+x2)
```
