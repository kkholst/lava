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
  label.width,
  contrast,
  null,
  level = NULL,
  type = NULL,
  var.adj = NULL,
  df = NULL,
  back.transform = NULL
)
```

## Arguments

- x:

  model object (`glm`, `lvmfit`, ...) or an existing `estimate` object.
  When two model objects are supplied (e.g., `estimate(g, g0)`) a
  likelihood-ratio test is performed.

- f:

  transformation of model parameters. Accepts several input types: - A
  **function** `f(p)` or `f(p, data)`: applies the delta method. When
  `f` returns a named list the names are used as parameter labels. - A
  **matrix**: used as a contrast (linear combination) matrix. - A
  **numeric vector** of parameter indices: converted to a contrast that
  selects and differences those parameters. - A **list** of indices:
  each element selects one parameter. - **Character** expressions:
  supports wildcards (`"?"`, `"*"`) and arithmetic on parameter names
  (e.g., `"z" - "x"`, `2 * "z" - 3 * "x"`).

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
  produce cluster-robust standard errors.

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
  of `f(p, data)` over all rows of `data`. The influence function
  accounts for both the empirical averaging and the parameter estimation
  uncertainty (see Details).

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

- contrast:

  (deprecated, use summary method) contrast matrix for a final Wald
  test. When supplied together with `null`, tests \\H_0: B\theta =
  b_0\\.

- null:

  (deprecated, use summary method) null hypothesis vector \\b_0\\ to
  test against (default 0)

- level:

  (deprecated, use summary method) level of confidence limits (default
  0.95)

- type:

  (deprecated, use summary method) type of small-sample correction for
  cluster-robust variance. One of: - `"robust"` (default): no
  correction.

  - `"df"`: applies \\n/(n-p)\\ correction (Mancl & DeRouen, 2001). -
    `"mbn"`: Morel-Bokossa-Neerchal (2003) correction. - `"hc3"`:
    leverage-adjusted HC3-type correction (blended with `var.adj`). -
    `"hc4"`: Cribari-Neto (2004) leverage-adjusted correction.

- var.adj:

  (deprecated, use summary method) blending parameter for the HC3
  leverage adjustment (default 0.25). Controls the weight between
  observation-level empirical leverage and the average leverage \\p/n\\.

- df:

  (deprecated, use summary method) degrees of freedom for t-based
  inference (default: `NULL` for Gaussian approximation; when set,
  confidence intervals and p-values use the t-distribution with `df`
  degrees of freedom).

- back.transform:

  (deprecated, use summary method) function applied to the point
  estimates and confidence interval bounds *after* inference is
  performed on the original scale. Useful for variance-stabilizing
  transformations, e.g., compute CIs on the `atanh` (Fisher z) scale and
  back-transform with `tanh`.

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
[pairwise.diff](https://kkholst.github.io/lava/reference/contr.md),
[c.estimate](https://kkholst.github.io/lava/reference/c.estimate.md),
[summary.estimate](https://kkholst.github.io/lava/reference/summary.estimate.md),
`coef.estimate`, `vcov.estimate`, `transform.estimate`,
`labels.estimate`,

## Examples

``` r
## Simulation from logistic regression model
m <- lvm(y~x+z);
distribution(m,y~x) <- binomial.lvm("logit")
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

## Marginalize / standardization
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

## Clusters and subset (conditional marginal effects)
d$id <- rep(seq(nrow(d)/4),each=4)
estimate(g,function(p,data)
         list(p0=expit(p[1] + p["z"]*data[,"z"])),
         subset=d$z>0, id=d$id, average=TRUE)
#>    Estimate Std.Err   2.5%  97.5%    P-value
#> p0   0.6754 0.02282 0.6307 0.7202 1.558e-192

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
#>   (Intercept)           x
#> 1 -0.04904705 -0.14189362
#> 1  0.14744213 -0.15724052
#> 4  1.15168804 -0.54934034
#> 1  0.37775868  0.06124426
#> 3 -0.01647188 -0.03981494
#> 1 -1.75658710  2.38806661
#> 2 -0.48187968 -0.18690182
#> 3 -0.56660619 -0.30858145
#> 4  0.77907752 -0.82047631
#> 5  0.41462553 -0.24506186
#> attr(,"bread")
#>             (Intercept)          x
#> (Intercept)   0.6136517 -0.1000796
#> x            -0.1000796  0.9011828
IC(estimate(l1,id=id1))
#>   (Intercept)           x
#> 1  -0.6402167  1.07508836
#> 4   0.9653828 -0.68490833
#> 3  -0.2915390 -0.17419820
#> 2  -0.2409398 -0.09345091
#> 5   0.2073128 -0.12253093

## Combined i.i.d. decomposition
e1 <- estimate(l1,id=id1)
e2 <- estimate(l2,id=id2)
e3 <- estimate(l3,id=id3)
(a2 <- merge(e1,e2,e3))
#>               Estimate Std.Err     2.5% 97.5%   P-value
#> (Intercept)     0.5418  0.2472  0.05731 1.026 0.0283949
#> x               0.9636  0.2592  0.45568 1.472 0.0002006
#> (Intercept).1   0.2156  0.5105 -0.78505 1.216 0.6728124
#> x.1             1.2369  0.4854  0.28554 2.188 0.0108277
#> (Intercept).2   0.4769  0.2979 -0.10705 1.061 0.1094554
#> x.2             1.1545  0.3714  0.42660 1.882 0.0018795

## If all models were estimated on the same data we could use the
## syntax:
## Reduce(merge,estimate(list(l1,l2,l3)))

## Same:
IC(a1 <- merge(l1,l2,l3,id=list(id1,id2,id3)))
#>    (Intercept)          x (Intercept).1         x.1 (Intercept).2         x.2
#> 1   -1.2804333  2.1501767    3.53262530 -3.15431642     0.7554987  2.18566521
#> 4    1.9307656 -1.3698167   -2.27428000 -0.36871845    -1.5561326 -0.25228854
#> 3   -0.5830781 -0.3483964    1.08830429 -0.51910711    -0.1215704  0.05798750
#> 2   -0.4818797 -0.1869018    0.06810688 -0.07263298     0.5200220 -0.55458052
#> 5    0.4146255 -0.2450619    0.25612485  0.61909106    -0.4880168 -1.17960761
#> 6    0.0000000  0.0000000   -2.62738927  3.57191544     1.6742259 -2.27609724
#> 7    0.0000000  0.0000000    0.33520853  0.13001396     0.8010558  0.31069745
#> 8    0.0000000  0.0000000   -0.37870058 -0.20624550     0.1460747  0.07955428
#> 9    0.0000000  0.0000000    0.00000000  0.00000000    -1.3102941  1.37992084
#> 10   0.0000000  0.0000000    0.00000000  0.00000000    -0.4208633  0.24874863

IC(merge(l1,l2,l3,id=TRUE)) # one-to-one (same clusters)
#>    (Intercept)           x (Intercept).1         x.1 (Intercept).2         x.2
#> 1  -0.04904705 -0.14189362   0.004626991  0.01338593     0.7554987  2.18566521
#> 2   0.14744213 -0.15724052   0.068106877 -0.07263298     0.5200220 -0.55458052
#> 3   1.15168804 -0.54934034   1.088304293 -0.51910711    -0.1215704  0.05798750
#> 4   0.37775868  0.06124426  -2.274280005 -0.36871845    -1.5561326 -0.25228854
#> 5  -0.01647188 -0.03981494   0.256124851  0.61909106    -0.4880168 -1.17960761
#> 6  -1.75658710  2.38806661  -2.627389265  3.57191544     1.6742259 -2.27609724
#> 7  -0.48187968 -0.18690182   0.335208530  0.13001396     0.8010558  0.31069745
#> 8  -0.56660619 -0.30858145  -0.378700581 -0.20624550     0.1460747  0.07955428
#> 9   0.77907752 -0.82047631   2.342596539 -2.46707796    -1.3102941  1.37992084
#> 10  0.41462553 -0.24506186   1.185401770 -0.70062440    -0.4208633  0.24874863
IC(merge(l1,l2,l3,id=FALSE)) # independence
#>    (Intercept)          x (Intercept).1         x.1 (Intercept).2        x.2
#> 1  -0.14714116 -0.4256809    0.00000000  0.00000000     0.0000000  0.0000000
#> 2   0.44232639 -0.4717216    0.00000000  0.00000000     0.0000000  0.0000000
#> 3   3.45506412 -1.6480210    0.00000000  0.00000000     0.0000000  0.0000000
#> 4   1.13327605  0.1837328    0.00000000  0.00000000     0.0000000  0.0000000
#> 5  -0.04941565 -0.1194448    0.00000000  0.00000000     0.0000000  0.0000000
#> 6  -5.26976131  7.1641998    0.00000000  0.00000000     0.0000000  0.0000000
#> 7  -1.44563904 -0.5607055    0.00000000  0.00000000     0.0000000  0.0000000
#> 8  -1.69981856 -0.9257444    0.00000000  0.00000000     0.0000000  0.0000000
#> 9   2.33723256 -2.4614289    0.00000000  0.00000000     0.0000000  0.0000000
#> 10  1.24387660 -0.7351856    0.00000000  0.00000000     0.0000000  0.0000000
#> 11  0.00000000  0.0000000    0.01388097  0.04015779     0.0000000  0.0000000
#> 12  0.00000000  0.0000000    0.20432063 -0.21789893     0.0000000  0.0000000
#> 13  0.00000000  0.0000000    3.26491288 -1.55732132     0.0000000  0.0000000
#> 14  0.00000000  0.0000000   -6.82284001 -1.10615534     0.0000000  0.0000000
#> 15  0.00000000  0.0000000    0.76837455  1.85727317     0.0000000  0.0000000
#> 16  0.00000000  0.0000000   -7.88216780 10.71574631     0.0000000  0.0000000
#> 17  0.00000000  0.0000000    1.00562559  0.39004188     0.0000000  0.0000000
#> 18  0.00000000  0.0000000   -1.13610174 -0.61873650     0.0000000  0.0000000
#> 19  0.00000000  0.0000000    7.02778962 -7.40123388     0.0000000  0.0000000
#> 20  0.00000000  0.0000000    3.55620531 -2.10187319     0.0000000  0.0000000
#> 21  0.00000000  0.0000000    0.00000000  0.00000000     2.2664960  6.5569956
#> 22  0.00000000  0.0000000    0.00000000  0.00000000     1.5600661 -1.6637416
#> 23  0.00000000  0.0000000    0.00000000  0.00000000    -0.3647111  0.1739625
#> 24  0.00000000  0.0000000    0.00000000  0.00000000    -4.6683977 -0.7568656
#> 25  0.00000000  0.0000000    0.00000000  0.00000000    -1.4640503 -3.5388228
#> 26  0.00000000  0.0000000    0.00000000  0.00000000     5.0226778 -6.8282917
#> 27  0.00000000  0.0000000    0.00000000  0.00000000     2.4031674  0.9320923
#> 28  0.00000000  0.0000000    0.00000000  0.00000000     0.4382241  0.2386628
#> 29  0.00000000  0.0000000    0.00000000  0.00000000    -3.9308824  4.1397625
#> 30  0.00000000  0.0000000    0.00000000  0.00000000    -1.2625898  0.7462459


# ------ influence function calculus -------
a <- estimate(coef = c("a" = 0.5), IC = scale(rnorm(10), scale=FALSE), id = 1:10)
b <- estimate(coef = c("b" = 0.8), IC = scale(rnorm(10), scale=FALSE), id = 1:10)

e <- c(a, b) # merge
merge(a, b)
#>   Estimate Std.Err   2.5%  97.5%  P-value
#> a      0.5  0.2346 0.0401 0.9599 0.033101
#> b      0.8  0.2868 0.2379 1.3621 0.005277
c(e1=a, b) # naming of par
#>    Estimate Std.Err   2.5%  97.5%  P-value
#> e1      0.5  0.2346 0.0401 0.9599 0.033101
#> b       0.8  0.2868 0.2379 1.3621 0.005277
labels(e, c("p1", "p2")) # renaming parameters
#>    Estimate Std.Err   2.5%  97.5%  P-value
#> p1      0.5  0.2346 0.0401 0.9599 0.033101
#> p2      0.8  0.2868 0.2379 1.3621 0.005277
e["a"] # subset
#>   Estimate Std.Err   2.5%  97.5% P-value
#> a      0.5  0.2346 0.0401 0.9599  0.0331
subset(e, "a")
#>   Estimate Std.Err   2.5%  97.5% P-value
#> a      0.5  0.2346 0.0401 0.9599  0.0331

# pipes
# c(a, b) |>
#  transform(function(x) x^2) |>
#  subset("a") |>
#  labels("sq")

# Parameter transformation with automatic calculation of derivatives
a * b
#>   Estimate Std.Err    2.5%  97.5% P-value
#> a      0.4  0.2583 -0.1063 0.9063  0.1215
(3 * cos(a) / sqrt(b) + 1) / a
#>   Estimate Std.Err   2.5% 97.5% P-value
#> a    7.887   4.783 -1.488 17.26 0.09917
expit(c(a,b))
#>   Estimate Std.Err   2.5%  97.5%   P-value
#> a   0.6225 0.05514 0.5144 0.7305 1.503e-29
#> b   0.6900 0.06134 0.5697 0.8102 2.382e-29
c(sum=sum(e), sum2=a+b,
  prod=prod(e), prod2=a*b)
#>       Estimate Std.Err    2.5%  97.5%  P-value
#> sum        1.3  0.4058  0.5047 2.0953 0.001356
#> sum2       1.3  0.4058  0.5047 2.0953 0.001356
#> prod       0.4  0.2583 -0.1063 0.9063 0.121526
#> prod2      0.4  0.2583 -0.1063 0.9063 0.121526
e %*% e # inner prod.
#>    Estimate Std.Err    2.5% 97.5% P-value
#> p1     0.89  0.5562 -0.2001  1.98  0.1096
c(1, 2) %*% e
#>    Estimate Std.Err   2.5% 97.5%  P-value
#> p1      2.1  0.6624 0.8018 3.398 0.001522
c(pow = a^b)
#>     Estimate Std.Err   2.5% 97.5%  P-value
#> pow   0.5743  0.2226 0.1382 1.011 0.009858
a^c(0.5, 2)
#>    Estimate Std.Err    2.5%  97.5%   P-value
#> p1   0.7071  0.1659  0.3819 1.0323 2.029e-05
#> p2   0.2500  0.2346 -0.2099 0.7099 2.867e-01
c(b=e["a"] * e["b"] / a, also.b=e["b"])
#>        Estimate Std.Err   2.5% 97.5%  P-value
#> b           0.8  0.2868 0.2379 1.362 0.005277
#> also.b      0.8  0.2868 0.2379 1.362 0.005277

B <- rbind(c(1,-1), c(1,0), c(0,1))
B %*% e
#>           Estimate Std.Err    2.5%  97.5%  P-value
#> [a] - [b]     -0.3  0.3316 -0.9499 0.3499 0.365625
#> a              0.5  0.2346  0.0401 0.9599 0.033101
#> b              0.8  0.2868  0.2379 1.3621 0.005277
e == 1 # wald-test, null-hypothesis H0: b=1
#> Call: estimate.default(data = NULL, id = id, coef = coefs, IC = ic0, 
#>     stack = FALSE, keep = keep)
#> ────────────────────────────────────────────────────────────
#>   Estimate Std.Err   2.5%  97.5% P-value
#> a      0.5  0.2346 0.0401 0.9599  0.0331
#> b      0.8  0.2868 0.2379 1.3621  0.4856
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [a] = 1
#>   [b] = 1 
#>  
#> chisq = 4.6135, df = 2, p-value = 0.09958
e == c(1,2)
#> Call: estimate.default(data = NULL, id = id, coef = coefs, IC = ic0, 
#>     stack = FALSE, keep = keep)
#> ────────────────────────────────────────────────────────────
#>   Estimate Std.Err   2.5%  97.5%   P-value
#> a      0.5  0.2346 0.0401 0.9599 3.310e-02
#> b      0.8  0.2868 0.2379 1.3621 2.859e-05
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [a] = 1
#>   [b] = 2 
#>  
#> chisq = 19.2204, df = 2, p-value = 6.704e-05
B %*% e == 1
#> Call: estimate.default(x = y, f = x)
#> ────────────────────────────────────────────────────────────
#>           Estimate Std.Err    2.5%  97.5%   P-value
#> [a] - [b]     -0.3  0.3316 -0.9499 0.3499 8.842e-05
#> a              0.5  0.2346  0.0401 0.9599 3.310e-02
#> b              0.8  0.2868  0.2379 1.3621 4.856e-01
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[a] - [b]] = 1
#>   [a] = 1
#>   [b] = 1 
#>  
#> chisq = 14.0808, df = 2, p-value = 0.0008758
```
