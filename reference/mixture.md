# Estimate mixture latent variable model.

Estimate mixture latent variable model

## Usage

``` r
mixture(
  x,
  data,
  k = length(x),
  control = list(),
  vcov = "observed",
  names = FALSE,
  ...
)
```

## Arguments

- x:

  List of `lvm` objects. If only a single `lvm` object is given, then a
  `k`-mixture of this model is fitted (free parameters varying between
  mixture components).

- data:

  `data.frame`

- k:

  Number of mixture components

- control:

  Optimization parameters (see details) \#type Type of EM algorithm
  (standard, classification, stochastic)

- vcov:

  of asymptotic covariance matrix (NULL to omit)

- names:

  If TRUE returns the names of the parameters (for defining starting
  values)

- ...:

  Additional arguments parsed to lower-level functions

## Details

Estimate parameters in a mixture of latent variable models via the EM
algorithm.

The performance of the EM algorithm can be tuned via the `control`
argument, a list where a subset of the following members can be altered:

- start:

  Optional starting values

- nstart:

  Evaluate `nstart` different starting values and run the EM-algorithm
  on the parameters with largest likelihood

- tol:

  Convergence tolerance of the EM-algorithm. The algorithm is stopped
  when the absolute change in likelihood and parameter (2-norm) between
  successive iterations is less than `tol`

- iter.max:

  Maximum number of iterations of the EM-algorithm

- gamma:

  Scale-down (i.e. number between 0 and 1) of the step-size of the
  Newton-Raphson algorithm in the M-step

- trace:

  Trace information on the EM-algorithm is printed on every `trace`th
  iteration

Note that the algorithm can be aborted any time (C-c) and still be saved
(via on.exit call).

## See also

`mvnmix`

## Author

Klaus K. Holst

## Examples

``` r

if (FALSE)  ## Reduce Ex.timings
m0 <- lvm(list(y~x+z,x~z))
distribution(m0,~z) <- dist_bernoulli()
#> Error: object 'm0' not found
d <- sim(m0,2000,p=c("y~z"=2,"y~x"=1),seed=1)
#> Error: object 'm0' not found

## unmeasured confounder example
m <- baptize(lvm(y~x, x~1));
intercept(m,~x+y) <- NA

if (requireNamespace('mets', quietly=TRUE)) {
  set.seed(42)
  M <- mixture(m,k=2,data=d,control=list(trace=1,tol=1e-6))
  summary(M)
  lm(y~x,d)
  estimate(M,"y~x")
  ## True slope := 1
}
#> Error: object 'd' not found
 # \dontrun{}
```
