# Influence function based inference

Primary tool for obtaining parameter estimates with robust (sandwich)
standard errors, applying the delta method, and testing linear
hypotheses. The function returns an object of class `estimate` which
serves as a general container for parameter estimates and their
influence functions (IFs). Three calling conventions are supported:

## Usage

``` r
# Default S3 method
estimate(
  x = NULL,
  f = NULL,
  ...,
  data,
  id,
  coef,
  IC = TRUE,
  vcov,
  stack = TRUE,
  average = FALSE,
  subset,
  keep,
  use,
  regex = FALSE,
  ignore.case = FALSE,
  print = NULL,
  labels,
  label.width
)
```

## Arguments

- x:

  model object (`glm`, `lvmfit`, ...) or an existing `estimate` object.
  When two model objects are supplied (e.g., `estimate(g, g0)`) a
  likelihood-ratio test is performed.

- f:

  transformation of model parameters. Accepts several input types:

  - A **function** `f(p)` or `f(p, data)`: applies the delta method.
    When `f` returns a named list the names are used as parameter
    labels.

  - A **matrix**: used as a contrast (linear combination) matrix. - A
    **numeric vector** of parameter indices: converted to a contrast
    that selects and differences those parameters.

  - A **list** of indices: each element selects one parameter.

  - **Character** expressions: supports wildcards (`"?"`, `"*"`) and
    arithmetic on parameter names (e.g., `"z" - "x"`,
    `2 * "z" - 3 * "x"`).

- ...:

  additional arguments to lower level functions

- data:

  `data.frame` used by `f` when the transformation depends on covariates
  (see `average`). Defaults to `model.frame(x)`.

- id:

  (optional) cluster identifier. Can be a vector of cluster IDs, a
  one-sided formula (evaluated in `data`), a single character column
  name, or a logical scalar (`TRUE` for one-to-one matching, `FALSE` for
  independence). When supplied, the IF is aggregated within clusters to
  produce cluster-robust standard errors. When `average = TRUE`, `id`
  refers to the rows of `data` (default: `rownames(data)`), and the ids
  of an `estimate` object (`index(x)`) are used as is (non-overlapping
  ids are treated as independent observations), whereas the rows of
  other model objects are linked via the row names of `data`.
  `id = NULL` removes the id (index) and the row names of the influence
  function from the returned object.

- coef:

  (optional) named parameter vector. Used instead of `coef(x)` when
  constructing an `estimate` object without a model.

- IC:

  if `TRUE` (default) the influence function matrix is estimated and
  stored in the returned object (extract with the
  [IC](https://kkholst.github.io/lava/reference/IC.default.md) method).
  Can also be a user-supplied IF matrix (one row per observation, one
  column per parameter), which is used directly instead of estimating it
  from `x`.

- vcov:

  (optional) covariance matrix of parameter estimates, or a logical. If
  `TRUE`, [stats::vcov](https://rdrr.io/r/stats/vcov.html) is used to
  obtain the (model-based) covariance matrix from `x`, yielding
  non-robust standard errors. If a matrix is supplied it is used
  directly. When omitted or `FALSE`, robust standard errors are computed
  from the influence function.

- stack:

  if `TRUE` (default) the influence function contributions are summed
  within each cluster defined by `id`. Set to `FALSE` to keep the
  un-stacked (per-observation) decomposition.

- average:

  if `TRUE` the function computes the standardized (marginalized)
  estimate \\\hat\Psi = P_n f(X; \hat\theta)\\, i.e., the empirical mean
  of `f(p, data)`, as defined by the `f` argument, over all rows of
  `data`. The influence function accounts for both the empirical
  averaging and the parameter estimation uncertainty (see Details).

- subset:

  (optional) logical vector, expression evaluated in `data`, or column
  name. When used together with `average = TRUE`, the average is
  conditioned on the subpopulation where `subset` is `TRUE`, yielding a
  conditional marginalized estimate.

- keep:

  (optional) index of parameters to keep from final result. Accepts
  integer indices, character names, or (with `regex = TRUE`)
  perl-compatible regular expressions.

- use:

  (optional) index of parameters to use in calculations. The selected
  parameters are first extracted (via `keep`) and then the remaining
  arguments (`f`, `contrast`, etc.) are applied to this subset.

- regex:

  if `TRUE` use perl-compatible regular expressions for `keep` and `use`
  arguments

- ignore.case:

  ignore case in regular expressions

- print:

  (optional) custom print function for the resulting `estimate` object

- labels:

  (optional) character vector of coefficient names

- label.width:

  (optional) max display width of labels

## Value

Object of class `estimate` with the following elements:

- coef:

  Named vector of parameter estimates.

- vcov:

  Variance-covariance matrix.

- IC:

  Influence function matrix (observations x parameters).

- coefmat:

  Formatted coefficient table (estimate, std.err, confidence limits,
  p-value).

- id:

  Cluster/id variable used.

- ncluster:

  Number of clusters.

- n:

  Number of observations.

- compare:

  (When `null` or contrasts are specified) Wald test result.

## Details

- `estimate(x, ...)` – extract estimates from a model object

- `estimate(coef=, IC=, ...)` – construct from coefficients and IF
  matrix

- `estimate(coef=, vcov=, ...)` – construct from coefficients and
  covariance matrix

## Influence functions and robust standard errors

An estimator \\\widehat{\theta}\\ is *regular and asymptotically linear*
(RAL) when it admits the iid decomposition
\$\$\sqrt{n}(\widehat{\theta}-\theta) = \frac{1}{\sqrt{n}}\sum\_{i=1}^n
\mathrm{IC}(Z_i; P) + o_p(1)\$\$ where \\\mathrm{IC}\\ is the unique
*influence function* satisfying \\E\\\mathrm{IC}(Z; P)\\ = 0\\. By the
central limit theorem \$\$\sqrt{n}(\widehat{\theta}-\theta)
\overset{d}{\longrightarrow} N(0,\\ \mathrm{Var}\\\mathrm{IC}(Z;
P)\\)\$\$ and the asymptotic variance is consistently estimated by the
empirical variance of the plugin IF estimate, yielding robust (sandwich)
standard errors. The estimated IF can be extracted with the
[IC](https://kkholst.github.io/lava/reference/IC.default.md) method.

## Parameter transformations (delta method)

When `f` is a function \\\phi: R^p \to R^m\\, the delta method is
applied: \$\$\sqrt{n}\\\phi(\widehat{\theta}) - \phi(\theta)\\ =
\frac{1}{\sqrt{n}}\sum\_{i=1}^n \nabla\phi(\theta)\\\mathrm{IC}(Z_i;
P) + o_p(1)\$\$ Derivatives are computed numerically via
[numDeriv::jacobian](https://rdrr.io/pkg/numDeriv/man/jacobian.html)
unless the function returns an attribute `"grad"` with the analytic
Jacobian.

Alternatively, `estimate` objects support direct arithmetic operations
(e.g., `a * b`, `exp(a)`, `a^b`) which apply the delta method with
*exact* (analytical) derivatives computed automatically. This influence
function calculus allows building complex transformations from simple
building blocks without numerical differentiation. See the last example
section ("influence function calculus") and
[`vignette("influencefunction", package = "lava")`](https://kkholst.github.io/lava/articles/influencefunction.md)
for details.

## Averaging and marginalization

When `average = TRUE` and `f(p, data)` depends on covariates, the target
parameter is the standardized (marginalized) estimate \\\Psi =
E\\f(X;\theta)\\\\. The IF for the averaged estimate accounts for both
the empirical averaging and parameter estimation uncertainty:
\$\$\mathrm{IC}\_\Psi(Z; P) = f(X;\theta) - \Psi + \[E\nabla\_\theta
f(X;\theta)\]\\\phi(Z; P)\$\$ When `subset` is also specified, the
average is conditioned on the subpopulation, yielding a conditional
marginalized estimate.

The model may be estimated on a different (e.g., smaller) dataset than
`data`. The two terms of the IF are then aligned by id: each term is
zero for ids outside its own support and rescaled by the inverse
proportion of observed ids (as in
[merge.estimate](https://kkholst.github.io/lava/reference/merge.estimate.md)).
If there are no common ids the model estimate and `data` are treated as
independent.

## Cluster-robust standard errors

When `id` is supplied, the per-observation IF contributions are summed
within clusters (when `stack = TRUE`), producing the cluster-level IF
\\\widetilde{\mathrm{IC}}(Z_i; P) = \sum\_{k=1}^{N_i}
\frac{n}{N}\mathrm{IC}(Z\_{ik}; P)\\. The resulting variance estimate is
equivalent to the GEE working independence sandwich estimator.

For full theoretical background and worked examples see
[`vignette("influencefunction", package = "lava")`](https://kkholst.github.io/lava/articles/influencefunction.md).

## See also

[estimate.array](https://kkholst.github.io/lava/reference/estimate.array.md),
[merge.estimate](https://kkholst.github.io/lava/reference/merge.estimate.md),
[contr](https://kkholst.github.io/lava/reference/contr.md),
[parsedesign](https://kkholst.github.io/lava/reference/contr.md),
[pairwise_diff](https://kkholst.github.io/lava/reference/contr.md),
[c.estimate](https://kkholst.github.io/lava/reference/c.estimate.md),
[summary.estimate](https://kkholst.github.io/lava/reference/summary.estimate.md),
`coef.estimate`, `vcov.estimate`, `transform.estimate`,
`labels.estimate`,

## Examples

``` r

## Simulation from logistic regression model
m <- lvm(y~x+z);
distribution(m,y~x) <- dist_bernoulli("logit")
d <- sim(m,1000)
g <- glm(y~z+x,data=d,family=binomial())
g0 <- glm(y~1,data=d,family=binomial())

## LRT
estimate(g, g0)
#> 
#>  - Likelihood ratio test -
#> 
#> data:  
#> chisq = 209.91, df = 2, p-value < 2.2e-16
#> sample estimates:
#> log likelihood (model 1) log likelihood (model 2) 
#>                -567.6493                -672.6041 
#> 


estimate(g)
#>              Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept) -0.001888 0.09785 -0.1937 0.1899 9.846e-01
#> z            0.953974 0.08114  0.7949 1.1130 6.481e-32
#> x            1.009058 0.14720  0.7205 1.2976 7.135e-12

## Testing contrasts
summary(estimate(g), null=0)
#> Call: estimate.default(x = x)
#> ────────────────────────────────────────────────────────────
#>              Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept) -0.001888 0.09785 -0.1937 0.1899 9.846e-01
#> z            0.953974 0.08114  0.7949 1.1130 6.481e-32
#> x            1.009058 0.14720  0.7205 1.2976 7.135e-12
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [(Intercept)] = 0
#>   [z] = 0
#>   [x] = 0 
#>  
#> chisq = 180.732, df = 3, p-value < 2.2e-16
estimate(g, rbind(c(1,1,0), c(1,0,2)))
#>                      Estimate Std.Err   2.5% 97.5%   P-value
#> [(Intercept)] + [z]    0.9521  0.1260 0.7052 1.199 4.075e-14
#> [(Intercept)] + 2[x]   2.0162  0.2404 1.5451 2.487 4.939e-17
summary(estimate(g, rbind(c(1,1,0), c(1,0,2))), null=c(1,2))
#> Call: estimate.default(x = x, f = ..1)
#> ────────────────────────────────────────────────────────────
#>                      Estimate Std.Err   2.5% 97.5% P-value
#> [(Intercept)] + [z]    0.9521  0.1260 0.7052 1.199  0.7037
#> [(Intercept)] + 2[x]   2.0162  0.2404 1.5451 2.487  0.9462
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[(Intercept)] + [z]] = 1
#>   [[(Intercept)] + 2[x]] = 2 
#>  
#> chisq = 0.1447, df = 2, p-value = 0.9302
estimate(g, 2:3) ## same as cbind(0,1,-1)
#>           Estimate Std.Err    2.5%  97.5% P-value
#> [z] - [x] -0.05508  0.1537 -0.3563 0.2461    0.72
estimate(g, as.list(2:3)) ## same as rbind(c(0,1,0),c(0,0,1))
#>   Estimate Std.Err   2.5% 97.5%   P-value
#> z    0.954 0.08114 0.7949 1.113 6.481e-32
#> x    1.009 0.14720 0.7205 1.298 7.135e-12
## Alternative syntax
estimate(g, "z", "z"-"x", 2*"z"-3*"x")
#>             Estimate Std.Err    2.5%   97.5%   P-value
#> z            0.95397 0.08114  0.7949  1.1130 6.481e-32
#> [z] - [x]   -0.05508 0.15367 -0.3563  0.2461 7.200e-01
#> 2[z] - 3[x] -1.11922 0.43991 -1.9814 -0.2570 1.095e-02
estimate(g, "?")  ## Wildcards
#>           Estimate Std.Err    2.5%  97.5% P-value
#> [z] - [x] -0.05508  0.1537 -0.3563 0.2461    0.72
estimate(g, "*Int*", "z")
#>              Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept) -0.001888 0.09785 -0.1937 0.1899 9.846e-01
#> z            0.953974 0.08114  0.7949 1.1130 6.481e-32
summary(estimate(g, "1", "2"-"3"), null = c(0,1))
#> Call: estimate.default(x = x, f = "1", ..2)
#> ────────────────────────────────────────────────────────────
#>              Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept) -0.001888 0.09785 -0.1937 0.1899 9.846e-01
#> [z] - [x]   -0.055083 0.15367 -0.3563 0.2461 6.606e-12
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [(Intercept)] = 0
#>   [[z] - [x]] = 1 
#>  
#> chisq = 77.8686, df = 2, p-value < 2.2e-16
estimate(g, 2, 3)
#>   Estimate Std.Err   2.5% 97.5%   P-value
#> z    0.954 0.08114 0.7949 1.113 6.481e-32
#> x    1.009 0.14720 0.7205 1.298 7.135e-12

## Usual (non-robust) confidence intervals
estimate(g, vcov=TRUE)
#>              Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept) -0.001888 0.09817 -0.1943 0.1905 9.847e-01
#> z            0.953974 0.08318  0.7909 1.1170 1.892e-30
#> x            1.009058 0.14639  0.7221 1.2960 5.469e-12
estimate(g, vcov=vcov(g))
#>              Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept) -0.001888 0.09817 -0.1943 0.1905 9.847e-01
#> z            0.953974 0.08318  0.7909 1.1170 1.892e-30
#> x            1.009058 0.14639  0.7221 1.2960 5.469e-12

## Transformations
estimate(g, function(p) p[1]+p[2])
#>             Estimate Std.Err   2.5% 97.5%   P-value
#> (Intercept)   0.9521   0.126 0.7052 1.199 4.075e-14

## Multiple parameters
e <- estimate(g, function(p) c(p[1]+p[2], p[1]*p[2]))
e
#>                Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept)    0.952086 0.12596  0.7052 1.1990 4.075e-14
#> (Intercept).1 -0.001801 0.09335 -0.1848 0.1812 9.846e-01
vcov(e)
#>             (Intercept) (Intercept)
#> (Intercept)  0.01586612 0.008982930
#> (Intercept)  0.00898293 0.008714966

## Label new parameters
estimate(g, function(p) list("a1"=p[1]+p[2], "b1"=p[1]*p[2]))
#>     Estimate Std.Err    2.5%  97.5%   P-value
#> a1  0.952086 0.12596  0.7052 1.1990 4.075e-14
#> b1 -0.001801 0.09335 -0.1848 0.1812 9.846e-01
#'
## Multiple group
m <- lvm(y~x)
m <- baptize(m)
d2 <- d1 <- sim(m,50,seed=1)
e <- estimate(list(m,m),list(d1,d2))
estimate(e) ## Wrong
#>        Estimate Std.Err     2.5%  97.5%   P-value
#> y@1      0.1044 0.08277 -0.05785 0.2666 2.073e-01
#> y~x@1    0.9665 0.08727  0.79541 1.1375 1.677e-28
#> y~~y@1   0.6764 0.10629  0.46803 0.8847 1.977e-10
ee <- estimate(e, id=rep(seq(nrow(d1)), 2)) ## Clustered
ee
#>        Estimate Std.Err    2.5%  97.5%   P-value
#> y@1      0.1044  0.1171 -0.1250 0.3338 3.725e-01
#> y~x@1    0.9665  0.1234  0.7246 1.2084 4.859e-15
#> y~~y@1   0.6764  0.1503  0.3817 0.9710 6.814e-06
estimate(lm(y~x,d1))
#>             Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept)   0.1044  0.1171 -0.1251 0.3338 3.726e-01
#> x             0.9665  0.1234  0.7246 1.2084 4.853e-15

## Standardization (g-computation)
f <- function(p,data)
  list(p0=expit(p["(Intercept)"] + p["z"]*data[,"z"]),
       p1=expit(p["(Intercept)"] + p["x"] + p["z"]*data[,"z"]))
e <- estimate(g, f, average=TRUE)
e
#>    Estimate Std.Err   2.5%  97.5%    P-value
#> p0   0.5010 0.02140 0.4591 0.5429 3.025e-121
#> p1   0.7007 0.01973 0.6620 0.7393 3.516e-276
estimate(e,diff)
#>    Estimate Std.Err  2.5%  97.5%   P-value
#> p1   0.1997  0.0279 0.145 0.2543 8.194e-13
estimate(e,cbind(1,1))
#>             Estimate Std.Err  2.5% 97.5% P-value
#> [p0] + [p1]    1.202 0.03027 1.142 1.261       0

# g-computation on non-overlapping data:
d$id <- paste0("i", 1:nrow(d))
d$w <- rbinom(nrow(d), 1, 0.5)
d1 <- subset(d, w == 1)
g1 <- glm(y ~ x + z, data=d1, family=binomial)
e1 <- estimate(g1, id=d1$id)
estimate(g1, f, data=d, id="id", average=TRUE)
#>    Estimate Std.Err   2.5%  97.5%    P-value
#> p0   0.5329 0.03077 0.4726 0.5933  3.334e-67
#> p1   0.6715 0.02827 0.6161 0.7269 9.955e-125

## Clusters and subset (conditional marginal effects)
d$id <- rep(seq(nrow(d)/4),each=4)
estimate(g,function(p,data)
         list(p0=expit(p[1] + p["z"]*data[,"z"])),
         subset=d$z>0, id=d$id, average=TRUE)
#>    Estimate Std.Err   2.5%  97.5%    P-value
#> p0   0.6754 0.02282 0.6307 0.7202 1.558e-192

## Model estimated on a subset, standardized over the full data
g1 <- glm(y~z+x, data=subset(d, id<=100), family=binomial())
estimate(g1, function(p,data) expit(p[1] + p["z"]*data[,"z"]),
         data=d, id=d$id, average=TRUE)
#>     Estimate Std.Err   2.5% 97.5%   P-value
#> val   0.5034 0.03295 0.4389 0.568 1.045e-52

## More examples with clusters:
m <- lvm(c(y1,y2,y3)~u+x)
d <- sim(m,10)
l1 <- glm(y1~x,data=d)
l2 <- glm(y2~x,data=d)
l3 <- glm(y3~x,data=d)

## Some random id-numbers
id1 <- c(1,1,4,1,3,1,2,3,4,5)
id2 <- c(1,2,3,4,5,6,7,8,1,1)
id3 <- seq(10)

## Un-stacked and stacked i.i.d. decomposition
IC(estimate(l1,id=id1,stack=FALSE))
#>   (Intercept)             x
#> 1 -1.30176154  0.8248362274
#> 1 -2.76419480  1.2576262814
#> 4  3.35596412 -1.8626151766
#> 1  0.75706350  0.0853556992
#> 3 -0.01210074  0.0007053009
#> 1  0.01064696 -1.6201534419
#> 2  0.16791446  1.0037531050
#> 3 -0.48471522 -0.1789769666
#> 4  0.04808727  0.3358573490
#> 5  0.22309598  0.1536116222
#> attr(,"bread")
#>             (Intercept)          x
#> (Intercept)   2.1652951 -0.9286132
#> x            -0.9286132  0.9731302
IC(estimate(l1,id=id1))
#>   (Intercept)           x
#> 1 -1.64912294  0.27383238
#> 4  1.70202570 -0.76337891
#> 3 -0.24840798 -0.08913583
#> 2  0.08395723  0.50187655
#> 5  0.11154799  0.07680581
#> attr(,"bread")
#>             (Intercept)          x
#> (Intercept)   2.1652951 -0.9286132
#> x            -0.9286132  0.9731302
#> attr(,"N")
#> [1] 10

## Combined i.i.d. decomposition
e1 <- estimate(l1,id=id1)
e2 <- estimate(l2,id=id2)
e3 <- estimate(l3,id=id3)
(a2 <- merge(e1,e2,e3))
#>               Estimate Std.Err     2.5%  97.5%   P-value
#> (Intercept)     0.4583  0.4774 -0.47743 1.3939 3.371e-01
#> x               0.2966  0.1922 -0.08012 0.6733 1.228e-01
#> (Intercept).1   0.4713  0.3098 -0.13600 1.0786 1.282e-01
#> x.1             0.4436  0.1723  0.10586 0.7813 1.004e-02
#> (Intercept).2   1.3949  0.2434  0.91777 1.8720 1.004e-08
#> x.2             0.3491  0.2805 -0.20070 0.8989 2.133e-01

## If all models were estimated on the same data we could use the
## syntax:
## Reduce(merge,estimate(list(l1,l2,l3)))

## Same:
IC(a1 <- merge(l1,l2,l3,id=list(id1,id2,id3)))
#>    (Intercept)          x (Intercept).1         x.1 (Intercept).2         x.2
#> 1   -3.2982459  0.5476648   1.174570191  0.25785719   -1.04890318  0.66461738
#> 4    3.4040514 -1.5267578  -0.539603784 -0.06083804    0.28895722  0.03257870
#> 3   -0.4968160 -0.1782717  -2.030027495  1.12669858   -0.07702585  0.04275061
#> 2    0.1679145  1.0037531   1.433578876 -0.65223568    1.71291544 -0.77932549
#> 5    0.2230960  0.1536116  -0.956823688  0.05576920   -1.24251391  0.07242087
#> 6    0.0000000  0.0000000   0.006847534 -1.04199253    0.01108649 -1.68703664
#> 7    0.0000000  0.0000000  -0.003888183 -0.02324264    0.32350393  1.93383032
#> 8    0.0000000  0.0000000   0.915346549  0.33798392   -0.22438808 -0.08285339
#> 9    0.0000000  0.0000000   0.000000000  0.00000000   -0.05932593 -0.41435180
#> 10   0.0000000  0.0000000   0.000000000  0.00000000    0.31569386  0.21736943

IC(merge(l1,l2,l3,id=TRUE)) # one-to-one (same clusters)
#>    (Intercept)             x (Intercept).1         x.1 (Intercept).2
#> 1  -1.30176154  0.8248362274   0.307279834 -0.19470197   -1.04890318
#> 2  -2.76419480  1.2576262814   1.433578876 -0.65223568    1.71291544
#> 3   3.35596412 -1.8626151766  -2.030027495  1.12669858   -0.07702585
#> 4   0.75706350  0.0853556992  -0.539603784 -0.06083804    0.28895722
#> 5  -0.01210074  0.0007053009  -0.956823688  0.05576920   -1.24251391
#> 6   0.01064696 -1.6201534419   0.006847534 -1.04199253    0.01108649
#> 7   0.16791446  1.0037531050  -0.003888183 -0.02324264    0.32350393
#> 8  -0.48471522 -0.1789769666   0.915346549  0.33798392   -0.22438808
#> 9   0.04808727  0.3358573490  -0.022969224 -0.16042463   -0.05932593
#> 10  0.22309598  0.1536116222   0.890259581  0.61298379    0.31569386
#>            x.2
#> 1   0.66461738
#> 2  -0.77932549
#> 3   0.04275061
#> 4   0.03257870
#> 5   0.07242087
#> 6  -1.68703664
#> 7   1.93383032
#> 8  -0.08285339
#> 9  -0.41435180
#> 10  0.21736943
IC(merge(l1,l2,l3,id=FALSE)) # independence
#>    (Intercept)            x (Intercept).1         x.1 (Intercept).2        x.2
#> 1  -3.90528462  2.474508682    0.00000000  0.00000000    0.00000000  0.0000000
#> 2  -8.29258439  3.772878844    0.00000000  0.00000000    0.00000000  0.0000000
#> 3  10.06789236 -5.587845530    0.00000000  0.00000000    0.00000000  0.0000000
#> 4   2.27119050  0.256067098    0.00000000  0.00000000    0.00000000  0.0000000
#> 5  -0.03630222  0.002115903    0.00000000  0.00000000    0.00000000  0.0000000
#> 6   0.03194089 -4.860460326    0.00000000  0.00000000    0.00000000  0.0000000
#> 7   0.50374338  3.011259315    0.00000000  0.00000000    0.00000000  0.0000000
#> 8  -1.45414566 -0.536930900    0.00000000  0.00000000    0.00000000  0.0000000
#> 9   0.14426182  1.007572047    0.00000000  0.00000000    0.00000000  0.0000000
#> 10  0.66928794  0.460834867    0.00000000  0.00000000    0.00000000  0.0000000
#> 11  0.00000000  0.000000000    0.92183950 -0.58410592    0.00000000  0.0000000
#> 12  0.00000000  0.000000000    4.30073663 -1.95670704    0.00000000  0.0000000
#> 13  0.00000000  0.000000000   -6.09008249  3.38009575    0.00000000  0.0000000
#> 14  0.00000000  0.000000000   -1.61881135 -0.18251412    0.00000000  0.0000000
#> 15  0.00000000  0.000000000   -2.87047106  0.16730760    0.00000000  0.0000000
#> 16  0.00000000  0.000000000    0.02054260 -3.12597759    0.00000000  0.0000000
#> 17  0.00000000  0.000000000   -0.01166455 -0.06972793    0.00000000  0.0000000
#> 18  0.00000000  0.000000000    2.74603965  1.01395175    0.00000000  0.0000000
#> 19  0.00000000  0.000000000   -0.06890767 -0.48127388    0.00000000  0.0000000
#> 20  0.00000000  0.000000000    2.67077874  1.83895137    0.00000000  0.0000000
#> 21  0.00000000  0.000000000    0.00000000  0.00000000   -3.14670954  1.9938521
#> 22  0.00000000  0.000000000    0.00000000  0.00000000    5.13874631 -2.3379765
#> 23  0.00000000  0.000000000    0.00000000  0.00000000   -0.23107756  0.1282518
#> 24  0.00000000  0.000000000    0.00000000  0.00000000    0.86687167  0.0977361
#> 25  0.00000000  0.000000000    0.00000000  0.00000000   -3.72754173  0.2172626
#> 26  0.00000000  0.000000000    0.00000000  0.00000000    0.03325947 -5.0611099
#> 27  0.00000000  0.000000000    0.00000000  0.00000000    0.97051180  5.8014910
#> 28  0.00000000  0.000000000    0.00000000  0.00000000   -0.67316423 -0.2485602
#> 29  0.00000000  0.000000000    0.00000000  0.00000000   -0.17797778 -1.2430554
#> 30  0.00000000  0.000000000    0.00000000  0.00000000    0.94708158  0.6521083


# ------ influence function calculus -------
ic1 <- scale(rnorm(10), scale=FALSE)
a <- estimate(coef = c("a" = 0.5), IC = ic1, id = 1:10)
b <- estimate(coef = c("b" = 0.8), IC = ic1, id = 1:10)

e <- c(a, b) # merge
merge(a, b)
#>   Estimate Std.Err    2.5% 97.5%  P-value
#> a      0.5  0.2778 -0.0444 1.044 0.071841
#> b      0.8  0.2778  0.2556 1.344 0.003974
c(e1=a, b) # naming of par
#>    Estimate Std.Err    2.5% 97.5%  P-value
#> e1      0.5  0.2778 -0.0444 1.044 0.071841
#> b       0.8  0.2778  0.2556 1.344 0.003974
labels(e, c("p1", "p2")) # renaming parameters
#>    Estimate Std.Err    2.5% 97.5%  P-value
#> p1      0.5  0.2778 -0.0444 1.044 0.071841
#> p2      0.8  0.2778  0.2556 1.344 0.003974
e["a"] # subset
#>   Estimate Std.Err    2.5% 97.5% P-value
#> a      0.5  0.2778 -0.0444 1.044 0.07184
subset(e, "a")
#>   Estimate Std.Err    2.5% 97.5% P-value
#> a      0.5  0.2778 -0.0444 1.044 0.07184

# pipes
# c(a, b) |>
#  transform(function(x) x^2) |>
#  subset("a") |>
#  labels("sq")

# Parameter transformation with automatic calculation of derivatives
a * b
#>   Estimate Std.Err    2.5% 97.5% P-value
#> a      0.4  0.3611 -0.3077 1.108   0.268
(3 * cos(a) / sqrt(b) + 1) / a
#>   Estimate Std.Err   2.5% 97.5% P-value
#> a    7.887   6.297 -4.454 20.23  0.2104
expit(c(a,b))
#>   Estimate Std.Err   2.5%  97.5%   P-value
#> a   0.6225 0.06527 0.4945 0.7504 1.484e-21
#> b   0.6900 0.05942 0.5735 0.8064 3.550e-31
c(sum=sum(e), sum2=a+b,
  prod=prod(e), prod2=a*b)
#>       Estimate Std.Err    2.5% 97.5% P-value
#> sum        1.3  0.5555  0.2112 2.389 0.01928
#> sum2       1.3  0.5555  0.2112 2.389 0.01928
#> prod       0.4  0.3611 -0.3077 1.108 0.26796
#> prod2      0.4  0.3611 -0.3077 1.108 0.26796
e %*% e # inner prod.
#>    Estimate Std.Err    2.5% 97.5% P-value
#> p1     0.89  0.7222 -0.5254 2.305  0.2178
c(1, 2) %*% e
#>    Estimate Std.Err   2.5% 97.5% P-value
#> p1      2.1  0.8333 0.4668 3.733 0.01173
c(pow = a^b)
#>     Estimate Std.Err   2.5%  97.5%   P-value
#> pow   0.5743  0.1447 0.2908 0.8579 7.186e-05
a^c(0.5, 2)
#>    Estimate Std.Err    2.5%  97.5%   P-value
#> p1   0.7071  0.1964  0.3222 1.0921 0.0003179
#> p2   0.2500  0.2778 -0.2944 0.7944 0.3680875
c(b=e["a"] * e["b"] / a, also.b=e["b"])
#>        Estimate Std.Err   2.5% 97.5%  P-value
#> b           0.8  0.2778 0.2556 1.344 0.003974
#> also.b      0.8  0.2778 0.2556 1.344 0.003974

B <- rbind(c(1,-1), c(1,0), c(0,1))
B %*% e
#>           Estimate Std.Err    2.5%  97.5%  P-value
#> [a] - [b]     -0.3  0.0000 -0.3000 -0.300 0.000000
#> a              0.5  0.2778 -0.0444  1.044 0.071841
#> b              0.8  0.2778  0.2556  1.344 0.003974
e == 1 # wald-test, null-hypothesis H0: b=1
#> Call: estimate.default(data = NULL, id = id, coef = coefs, IC = ic0, 
#>     stack = FALSE, keep = keep)
#> ────────────────────────────────────────────────────────────
#>   Estimate Std.Err    2.5% 97.5% P-value
#> a      0.5  0.2778 -0.0444 1.044 0.07184
#> b      0.8  0.2778  0.2556 1.344 0.47149
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [a] = 1
#>   [b] = 1 
#>  
#> chisq = 1.5878, df = 1, p-value = 0.2076
e == c(1,2)
#> Call: estimate.default(data = NULL, id = id, coef = coefs, IC = ic0, 
#>     stack = FALSE, keep = keep)
#> ────────────────────────────────────────────────────────────
#>   Estimate Std.Err    2.5% 97.5%   P-value
#> a      0.5  0.2778 -0.0444 1.044 7.184e-02
#> b      0.8  0.2778  0.2556 1.344 1.558e-05
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [a] = 1
#>   [b] = 2 
#>  
#> chisq = 9.3649, df = 1, p-value = 0.002212
B %*% e == 1
#> Call: estimate.default(x = y, f = x)
#> ────────────────────────────────────────────────────────────
#>           Estimate Std.Err    2.5%  97.5% P-value
#> [a] - [b]     -0.3  0.0000 -0.3000 -0.300 0.00000
#> a              0.5  0.2778 -0.0444  1.044 0.07184
#> b              0.8  0.2778  0.2556  1.344 0.47149
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[a] - [b]] = 1
#>   [a] = 1
#>   [b] = 1 
#>  
#> chisq = 1.5878, df = 1, p-value = 0.2076
```
