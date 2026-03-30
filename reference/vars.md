# Extract variable names from latent variable model

Extract exogenous variables (predictors), endogenous variables
(outcomes), latent variables (random effects), manifest (observed)
variables from a `lvm` object.

## Usage

``` r
vars(x,...)

endogenous(x,...)

exogenous(x,...)

manifest(x,...)

latent(x,...)

# S3 method for class 'lvm'
exogenous(x, xfree = TRUE, ...) <- value

# S3 method for class 'lvm'
exogenous(x,variable,latent=FALSE,index=TRUE,...)

# S3 method for class 'lvm'
latent(x, clear = FALSE, ...) <- value
```

## Arguments

- x:

  `lvm`-object

- ...:

  Additional arguments to be passed to the low level functions

- variable:

  list of variables to alter

- latent:

  Logical defining whether latent variables without parents should be
  included in the result

- index:

  For internal use only

- clear:

  Logical indicating whether to add or remove latent variable status

- xfree:

  For internal use only

- value:

  Formula or character vector of variable names.

## Value

Vector of variable names.

## Details

`vars` returns all variables of the `lvm`-object including manifest and
latent variables. Similarily `manifest` and `latent` returns the
observered resp. latent variables of the model. `exogenous` returns all
manifest variables without parents, e.g. covariates in the model,
however the argument `latent=TRUE` can be used to also include latent
variables without parents in the result. Pr. default `lava` will not
include the parameters of the exogenous variables in the optimisation
routine during estimation (likelihood of the remaining observered
variables conditional on the covariates), however this behaviour can be
altered via the assignment function `exogenous<-` telling `lava` which
subset of (valid) variables to condition on. Finally `latent` returns a
vector with the names of the latent variables in `x`. The assigment
function `latent<-` can be used to change the latent status of variables
in the model.

## See also

`endogenous`, `manifest`, `latent`, `exogenous`, `vars`

## Author

Klaus K. Holst

## Examples

``` r
g <- lvm(eta1 ~ x1+x2)
regression(g) <- c(y1,y2,y3) ~ eta1
latent(g) <- ~eta1
endogenous(g)
#> [1] "y1" "y2" "y3"
exogenous(g)
#> [1] "x1" "x2"
identical(latent(g), setdiff(vars(g),manifest(g)))
#> [1] TRUE
```
