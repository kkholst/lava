# Add covariance structure to Latent Variable Model

Define covariances between residual terms in a `lvm`-object.

## Usage

``` r
# S3 method for class 'lvm'
covariance(object, var1 = NULL, var2 = NULL, constrain = FALSE, pairwise = FALSE, ...) <- value
```

## Arguments

- object:

  `lvm`-object

- ...:

  Additional arguments to be passed to the low level functions

- var1:

  Vector of variables names (or formula)

- var2:

  Vector of variables names (or formula) defining pairwise covariance
  between `var1` and `var2`)

- constrain:

  Define non-linear parameter constraints to ensure positive definite
  structure

- pairwise:

  If TRUE and `var2` is omitted then pairwise correlation is added
  between all variables in `var1`

- value:

  List of parameter values or (if `var1` is unspecified)

## Value

A `lvm`-object

## Details

The `covariance` function is used to specify correlation structure
between residual terms of a latent variable model, using a formula
syntax.

For instance, a multivariate model with three response variables,

\$\$Y_1 = \mu_1 + \epsilon_1\$\$

\$\$Y_2 = \mu_2 + \epsilon_2\$\$

\$\$Y_3 = \mu_3 + \epsilon_3\$\$

can be specified as

`m <- lvm(~y1+y2+y3)`

Pr. default the two variables are assumed to be independent. To add a
covariance parameter \\r = cov(\epsilon_1,\epsilon_2)\\, we execute the
following code

`covariance(m) <- y1 ~ f(y2,r)`

The special function `f` and its second argument could be omitted thus
assigning an unique parameter the covariance between `y1` and `y2`.

Similarily the marginal variance of the two response variables can be
fixed to be identical (\\var(Y_i)=v\\) via

`covariance(m) <- c(y1,y2,y3) ~ f(v)`

To specify a completely unstructured covariance structure, we can call

`covariance(m) <- ~y1+y2+y3`

All the parameter values of the linear constraints can be given as the
right handside expression of the assigment function `covariance<-` if
the first (and possibly second) argument is defined as well. E.g:

`covariance(m,y1~y1+y2) <- list("a1","b1")`

`covariance(m,~y2+y3) <- list("a2",2)`

Defines

\$\$var(\epsilon_1) = a1\$\$

\$\$var(\epsilon_2) = a2\$\$

\$\$var(\epsilon_3) = 2\$\$

\$\$cov(\epsilon_1,\epsilon_2) = b1\$\$

Parameter constraints can be cleared by fixing the relevant parameters
to `NA` (see also the `regression` method).

The function `covariance` (called without additional arguments) can be
used to inspect the covariance constraints of a `lvm`-object.

## See also

`regression<-`, `intercept<-`, `constrain<-` `parameter<-`, `latent<-`,
`cancel<-`, `kill<-`

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm()
### Define covariance between residuals terms of y1 and y2
covariance(m) <- y1~y2
covariance(m) <- c(y1,y2)~f(v) ## Same marginal variance
covariance(m) ## Examine covariance structure
#> Covariance parameters:
#>       y1 y2
#>    y1 v  * 
#>    y2 *  v 

```
