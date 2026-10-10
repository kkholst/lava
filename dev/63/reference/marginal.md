# Marginal (standardized) estimates via g-computation

Estimates the marginal (standardized) parameter \\\Psi = E\\f(X;
\theta)\\\\ by the empirical average of `f(p, data)` over the rows of
`data` (g-computation), where the parameter \\\theta\\ is estimated by
`object`. The influence function of the estimate accounts for both the
empirical averaging and the uncertainty of \\\widehat\theta\\.

## Usage

``` r
marginal(object, f, data, id, subset, labels = NULL, ...)
```

## Arguments

- object:

  model object (`glm`, `lvmfit`, ...) or `estimate` object providing the
  parameter estimates and their influence function.

- f:

  function where the first argument is the parameter vector and the
  optional named arguments `data` and `object` (the model object, e.g.
  an `estimate` object. Should return a vector, matrix or a (named) list
  of vectors with one value per observation (row of `data`). If the
  result has an attribute `"grad"` it is used as the Jacobian of `f`
  with respect to `p`

- data:

  `data.frame` over which `f` is averaged. Defaults to
  `model.frame(object)` (which is not available for `estimate` objects).

- id:

  (optional) ids of the values returned by `f` (the g-computation part):
  a vector with one value per value of `f`, or a column name or
  one-sided formula evaluated in `data` (requiring one value of `f` per
  row of `data`). Defaults to the names of the values of `f`, and
  otherwise `rownames(data)`. `id = NULL` uses the default ids but
  removes the id (index) from the returned object.

- subset:

  (optional) logical vector (one value per value of `f`), expression
  evaluated in `data` (columns of `data` take precedence over variables
  in the calling environment), or column name. The average is then
  conditioned on the subpopulation where `subset` is `TRUE` (conditional
  marginal estimate).

- labels:

  (optional) character vector of parameter names

- ...:

  additional arguments passed to `f`

## Value

`estimate` object

## Details

The influence function of the estimate is \$\$\mathrm{IC}\_\Psi(Z; P) =
f(X;\theta) - \Psi + \[E\nabla\_\theta f(X;\theta)\]\\\phi(Z; P)\$\$
where \\\phi\\ is the influence function of \\\widehat\theta\\.

The two terms are identified by different ids: the first term by `id`
(the rows of `data`), and the second term by the ids of the model, i.e.,
`index(object)` for `estimate` objects, or the row names of the
influence function (for model objects such as `glm` the row names of the
model frame). The model may therefore be estimated on another (e.g.,
smaller or partly overlapping) dataset than `data`. The two terms are
aligned on the union of ids: each term is zero for ids outside its own
support and rescaled by the inverse proportion of observed ids (as in
[merge.estimate](https://kkholst.github.io/lava/reference/merge.estimate.md)).
If there are no common ids, the model estimate and `data` are treated as
independent. Use `estimate(object, id=...)` to assign ids (or clusters)
to the model.

## See also

[estimate.default](https://kkholst.github.io/lava/reference/estimate.default.md),
[`vignette("influencefunction", package = "lava")`](https://kkholst.github.io/lava/articles/influencefunction.md)

## Examples

``` r
m <- lvm(y ~ x + z)
distribution(m, y ~ x) <- dist_bernoulli("logit")
d <- sim(m, 1000, seed = 1)
g <- glm(y ~ z + x, data = d, family = binomial())

## Standardization (g-computation)
f <- function(p, data)
  list(p0 = expit(p["(Intercept)"] + p["z"] * data[, "z"]),
       p1 = expit(p["(Intercept)"] + p["x"] + p["z"] * data[, "z"]))
e <- marginal(g, f)
e
#>    Estimate Std.Err   2.5%  97.5%    P-value
#> p0   0.5143 0.02111 0.4729 0.5557 4.099e-131
#> p1   0.7097 0.02003 0.6705 0.7490 4.446e-275
estimate(e, diff)
#>    Estimate Std.Err   2.5% 97.5%   P-value
#> p1   0.1954 0.02788 0.1408  0.25 2.402e-12

## Model estimated on a subset, standardized over the full data
d$id <- paste0("i", seq_len(nrow(d)))
d$w <- rbinom(nrow(d), 1, 0.5)
d1 <- subset(d, w == 1)
g1 <- glm(y ~ x + z, data = d1, family = binomial)
e1 <- estimate(g1, id = d1$id)
marginal(e1, f, data = d, id = "id")
#>    Estimate Std.Err   2.5%  97.5%    P-value
#> p0   0.5118 0.02939 0.4542 0.5694  6.388e-68
#> p1   0.6927 0.02776 0.6383 0.7471 2.069e-137

## Conditional marginal effects (subset) and clusters
d$cl <- rep(seq(nrow(d) / 4), each = 4)
marginal(estimate(g, id = d$cl),
         function(p, data) expit(p[1] + p["z"] * data[, "z"]),
         data = d, subset = z > 0, id = "cl")
#>    Estimate Std.Err   2.5%  97.5%    P-value
#> p1   0.6903 0.02496 0.6414 0.7392 2.202e-168
```
