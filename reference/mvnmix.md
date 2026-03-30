# Estimate mixture latent variable model

Estimate mixture latent variable model

## Usage

``` r
mvnmix(
  data,
  k = 2,
  theta,
  steps = 500,
  tol = 1e-16,
  lambda = 0,
  mu = NULL,
  silent = TRUE,
  extra = FALSE,
  n.start = 1,
  init = "kmpp",
  ...
)
```

## Arguments

- data:

  `data.frame`

- k:

  Number of mixture components

- theta:

  Optional starting values

- steps:

  Maximum number of iterations

- tol:

  Convergence tolerance of EM algorithm

- lambda:

  Regularisation parameter. Added to diagonal of covariance matrix (to
  avoid singularities)

- mu:

  Initial centres (if unspecified random centres will be chosen)

- silent:

  Turn on/off output messages

- extra:

  Extra debug information

- n.start:

  Number of restarts

- init:

  Function to choose initial centres

- ...:

  Additional arguments parsed to lower-level functions

## Value

A `mixture` object

## Details

Estimate parameters in a mixture of latent variable models via the EM
algorithm.

## See also

`mixture`

## Author

Klaus K. Holst

## Examples

``` r
data(faithful)
set.seed(1)
M1 <- mvnmix(faithful[,"waiting",drop=FALSE],k=2)
M2 <- mvnmix(faithful,k=2)
if (interactive()) {
    par(mfrow=c(2,1))
    plot(M1,col=c("orange","blue"),ylim=c(0,0.05))
    plot(M2,col=c("orange","blue"))
}
```
