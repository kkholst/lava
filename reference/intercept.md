# Fix mean parameters in 'lvm'-object

Define linear constraints on intercept parameters in a `lvm`-object.

## Usage

``` r
# S3 method for class 'lvm'
intercept(object, vars, ...) <- value
```

## Arguments

- object:

  `lvm`-object

- ...:

  Additional arguments

- vars:

  character vector of variable names

- value:

  Vector (or list) of parameter values or labels (numeric or character)
  or a formula defining the linear constraints (see also the
  `regression` or `covariance` methods).

## Value

A `lvm`-object

## Details

The `intercept` function is used to specify linear constraints on the
intercept parameters of a latent variable model. As an example we look
at the multivariate regression model

\$\$ E(Y_1\|X) = \alpha_1 + \beta_1 X\$\$ \$\$ E(Y_2\|X) = \alpha_2 +
\beta_2 X\$\$

defined by the call

`m <- lvm(c(y1,y2) ~ x)`

To fix \\\alpha_1=\alpha_2\\ we call

`intercept(m) <- c(y1,y2) ~ f(mu)`

Fixed parameters can be reset by fixing them to `NA`. For instance to
free the parameter restriction of \\Y_1\\ and at the same time fixing
\\\alpha_2=2\\, we call

`intercept(m, ~y1+y2) <- list(NA,2)`

Calling `intercept` with no additional arguments will return the current
intercept restrictions of the `lvm`-object.

## Note

Variables will be added to the model if not already present.

## See also

`covariance<-`, `regression<-`, `constrain<-`, `parameter<-`,
`latent<-`, `cancel<-`, `kill<-`

## Author

Klaus K. Holst

## Examples

``` r

## A multivariate model
m <- lvm(c(y1,y2) ~ f(x1,beta)+x2)
regression(m) <- y3 ~ f(x1,beta)
intercept(m) <- y1 ~ f(mu)
intercept(m, ~y2+y3) <- list(2,"mu")
intercept(m) ## Examine intercepts of model (NA translates to free/unique paramete##r)
#> Intercept parameters:
#>     y1 y2 y3
#>     *  2  mu

```
