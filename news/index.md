# Changelog

## lava (development version)

- `merge.estimate` by default no longer sorts IC by id (old behaviour
  via new argument `sort=TRUE`)

## lava 1.9.2

CRAN release: 2026-06-30

- `estimate.default`: the `null`, `contrast`, `type`, and `var.adj`
  arguments are soft-deprecated. Use
  `summary(estimate(...), null=, contrast=, type=, var.adj=)` instead.
  The default-path Wald p-value (H0: beta = 0) continues to be reported
  by
  [`estimate()`](https://kkholst.github.io/lava/reference/estimate.default.md).
- `estimate.default`: removed `R`, `null.sim`, `score.deriv` and `folds`
  arguments.
- new `c.summary.estimate` S3-method for concatenating
  `summary.estimate` objects
- sim.default can now operate on function return objects constructed via
  `c.estimate(estimate_object, extra_args)` or
  `c.estimate(summary(estimate_object), extra_args)`
- `predict_glm` function
- fixed issue when package formula.tools was loaded due to overwriting
  of as.character.formula.
- bug-fix: `regression(model, "y", "x", value="b")` now works as
  expected
- `index.estimate`, `index<-.estimate` methods for getting and setting
  id/cluster

## lava 1.9.1

CRAN release: 2026-05-14

- Safe evaluation of rank in `wald_test`
- adding CI Length to `summary.sim` output
- fixing bug related to `estimate.index` in `summary.sim`
- updated `plot.sim`, `forestplot`, `plot.estimate`
- `estimate.default`: small-sample `type="hc3"` variance estimates
- `merge`, `c.estimate`: new `drop.ic` argument for dropping influence
  functions before merging
- `sim.default`: export seeds to replicate results
- bug-fix: `estimate.lvm` fixed issue with interval censored
  observations
- deprecated: `only.coef` argument in
  [`estimate.default()`](https://kkholst.github.io/lava/reference/estimate.default.md).
  Use `parameter(estimate(...))` instead.

## lava 1.9.0

CRAN release: 2026-04-05

- `estimate`: estimate objects can now be transformed via functions:
  `log`, `exp`, `+`, `-`, `*`, `/`, `^`, … See the influence vignette
  for full details and list of available mathematical transformations.
- `sim.default`: the simulation function `f`’s return object can now be
  an `estimate` object
- `summary.sim`: automatically derives estimate, se, confint parameters
  if the sim routine returns an `estimate` object
- `merge.estimate`: cast warning if `back.transform` was used
- `merge.estimate`: works with objects with and without influence
  function
- `summary.sim` `df` argument for calculating CIs based on t-dist.
  approximation
- breaking change: `+` operator for estimate objects no longer merges
  objects, use `%++%` or `c(...)` instead
- breaking change: `addattr`, `intfix`, `covfix` no longer exported.

## lava 1.8.2

CRAN release: 2025-10-30

- Improved closed testing procedure `closed_testing` (depr.
  `closed.testing`)
- More tests and documentation
- New data `deprdiag`

## lava 1.8.1

CRAN release: 2025-01-12

- `sim.default` now accepts the argument `R` to be a list (of lists) of
  arguments.
- New methods `subset.estimate`, `transform.estimate`, `labels.estimate`

## lava 1.8.0

CRAN release: 2024-03-05

- New methods `estimate.mlm`, `IC.mlm`, `pars.mlm`, `estimate.array`,
  `estimate.data.frame`
- `Print` method for tabular data (matrix, data.frame, data.table)
- `merge` now supports regular expressions
- `IC` returns row-names (default id) as obtained from model.matrix or
  similar
- New vignette: Influence Functions
- operators `%in.open%`, `%in.closed%` for checking if elements are
  within a range `3 %in.open% c(0,1)`)
- `estimate(..., estimator='glm')` now works with formulas with just an
  intercept
- `as.data.frame.sim`, `as.matrix.sim`
- fixed issues with quasi\* families and negative binomial regression
  (`MASS:glm.nb`)

## lava 1.7.3

CRAN release: 2023-11-04

- `parameter.estimate` method to extract matrix with estimates, standard
  errors, and confidence limits from and estimate object (coefmat
  element)
- pairwise difference with `'-'.estimate` and `pairwise.diff`
- optional mc.cores arguments to `cv` and `bootstrap`
- `parameter.lvm` now automatically removes previously variables in the
  lvm object with same name as new added parameters.
- `Print` function deals more gracefully with non-rectangular objects
- bug-fix in `stack.estimate` (wrong stand-errors in `twostage` since
  version 1.7.0)

## lava 1.7.2.1

CRAN release: 2023-02-27

- Maintenance release as version 1.7.2 broke compatibility with R\<4.1.

## lava 1.7.2

CRAN release: 2023-02-23

- Compatibility issues with development version of R fixed.
- cluster.index now also works when not loading the package (directly
  calling lava::estimate)
- `weibull.lvm` and `coxExponential.lvm` now uses default
  parametrizations similar to `rweibull`, `rexp`. `weibull.lvm` now has
  arguments “intercept”,“sigma” that directly relates to the accelerated
  failure time formulation.
- Packages `gof`, `lava.tobit` are removed from Suggested packages.

## lava 1.7.1

CRAN release: 2023-01-06

- Fixed bug in variance estimates from `estimate` with clustered
  observations.
- Discrete uniform distributions can now be specified with
  `uniform.lvm(value=...)`.

## lava 1.7.0

CRAN release: 2022-10-25

- `cv` method moved to the ‘targeted’ package
- New `IC` method that returns influence function of a model object. The
  `iid` argument `iid` of the `estimate` method is now replaced with an
  argument `IC` (with a user supplied matrix this must now be the actual
  influence function and not the sample-size scaled version returned by
  the `iid` method).
- fixed bug where calls like `regression("y", value=function(x) x)` did
  not work. `merge.estimate` now works without IC element

## lava 1.6.10

CRAN release: 2021-09-02

- Improved starting values for MLE optimization.
- New simulation distributions: `multinomial.lvm`, `none.lvm`,
  `constant.lvm`, `id.lvm`.
- `regression`, `regression.lvm`: the ‘value’ argument can now be a
  (non-linear) function specifying the functional relationship between
  outcomes and covariates (for simulation with the `sim` method).
- New `intervention` method for applying interventions on `lvm`-objects
- Progress updates are now done via the `progressr` library (enabled
  with `progressr::handlers(global=TRUE)`).
- Parallelization is now controlled via the future library. To enable
  multicore parallelization: `future::plan("multicore")`.
- New `plot_region` function for adding confidence regions to plots.

## lava 1.6.9

CRAN release: 2021-03-11

- `idplot`: now accepts matrix or data.frame as 1st argument. New
  argument: return.data.
- Unit tests updated
- Bug fixes: `cv`: rmse output fixed. score: Fixed bug for linear
  Gaussian model with argument ‘indiv=TRUE’. estimate.formula: call
  object initialized correctly. `plot.lvm`: ‘noplot’ argument now works
  with all plot engines.

## lava 1.6.8.1

CRAN release: 2020-11-04

- Maintenance release
- `confpred`: split-conformal prediction method updated

## lava 1.6.8

CRAN release: 2020-09-26

- Bug-fix: `parameter(m,x)` now returns an lvm object and not just x
- profile likelihood confidence intervals with tobit/censored
  observations.
- Vignettes added
  - Estimating partial correlations
  - Non-linear latent variable omdels
- Pseudo-inverse used with “normal” estimator
- Starting values for mixture fixed

## lava 1.6.7

CRAN release: 2020-03-05

- Fixed bug in the composite likelihood `complik` when used with
  censored variables (Surv objects).
- Fixed regular expression in ‘spaghetti’ function
- `plot.sim`: ‘rug’ argument is now by default FALSE and ‘auto.layout’
  disabled when nr=nc=1.
- base::sequence() is now a generic function and as a consequence
  `sequence.lvm` has been renamed to `Sequence.lvm`. The function
  `binary.lvm` is now an alias of `ones.lvm`.

## lava 1.6.6

CRAN release: 2019-08-01

- Weighted kmeans++ (`wkm`). Gaussian mixture models (`mvnmix`) are now
  initialized by default using kmeans++.
- `sim` method implemented for mvnmix models.
- Bug fix: Newton-Raphson method
  ([`lava::NR`](https://kkholst.github.io/lava/reference/NR.md)) used a
  numerical approximation of the Hessian even when submitted as
  attribute to the objective function.

## lava 1.6.5

CRAN release: 2019-02-12

- Maintenance release.

## lava 1.6.4

CRAN release: 2018-11-25

- New simulation distributions: constant relative risk and risk
  difference models as in Richardson, Robins and Wang, 2017):
  `binomial.rd`, `binomial.rr`. Base on new hook
  ‘simulate_multiple_inputs’ which allows the distribution to depend
  non-linearly on multiple different input variables.
- `sim.lvm`: ‘X’ argument can now fix (manipulate) any variable and not
  only exogenous variables.
- Summary function for `sim.default` updated (the ‘estimate’ argument
  can now be a list with each element being the estimate position and
  optionally standard error and true value).
- Starting values updated for mixture models. The parameter names can be
  obtained with `mixture(...,names=TRUE)` and set with
  `mixture(...,control=list(start=...)))`.
- Naming conventions for multigroup parameters: ‘<par@g>’ (par: name of
  parameter, g: first group number where ‘par’ is observed). Starting
  values can be specified with estimate(…,control(list(start=…))).
- New print and summary methods for mixture models.
- Renamed (weighted) K-means function ‘km’ to `wkm`.
- Derivative method deriv.function based on complex step derivatives.
- `twostageCV`: estimation of mixture models are now parallelized if
  mc.cores\>1.

## lava 1.6.3

CRAN release: 2018-08-10

- Fixed problems with plots (Rgraphviz)
- Better print method for `twostageCV`
- Improved M-step in `mixture` method

## lava 1.6.2

CRAN release: 2018-07-02

- `twostageCV`: cross-validate two-stage estimator
- rmvn, dmvn moved to mets package (C++ implementation, old versions
  renamed to
  [`lava::rmvn0`](https://kkholst.github.io/lava/reference/internal.md),
  [`lava::dmvn0`](https://kkholst.github.io/lava/reference/internal.md))
- mediation proportion handled correctly when direct effect is zero
- unit tests clean-up (namespace)
- `merge.lvm` now correctly handles fixed covariance parameters

## lava 1.6.1

CRAN release: 2018-03-28

- Newton-raphson algorithm made more robust.
- New `sim.as` method. `plot.sim` method now by default only plots
  density estimates
- Compatibility fix with Matrix library

## lava 1.6

CRAN release: 2018-01-12

- Mixture Latent variable models (`mixture`). Fast version requires
  ‘mets’ packages; Gaussian mixture models (`mvnmix`); weighted k-means
  (`km`)
- `estimate.default`: ‘keep’, ‘use’ arguments can be specified as
  regular expressions (with argument regex=TRUE). Summary method now
  returns Wald test (null: all parameters being zero).
- `makemissing`: seed argument added.
- Global change: ‘silent’ argument renamed to ‘messages’
- New utility functions: `Grep`, `Na2x`, `x2NA`, `wait`, `waitclick`,
  `rotation`, `Rot2d`, `Rot3d`
- Condition numbers calculated via SVD
- `na.pass0`: returns data.frame with original number of rows but with
  zeros (or first level of factors) in the rows with missing data.
- `stack`: ‘weights’ argument renamed to ‘propensity’. If
  propensity=TRUE, the first argument (model) will be treated as
  propensity score model (glm) and ‘predict’ method will be used for the
  predictions.
- `estimate.formula` now by default wraps glm such that the `iid` method
  return matrix of same size as full data (with zero rows where data are
  missing).
- Updated output functions for class ‘sim’ (print method and plot)..
  Plot method: density.alpha is applied to each standard error (‘se’)
  level.
- composite likelihood (`complik`) refactored + new example. `ordinal`
  method now cleans up properly when variables are removed (`rmvar`,
  `subset`).
- `twostage`: fixed for mixture model (class ‘lvm.mixture’). New help
  page + examples. Predict function updated (newdata argument where
  covariate levels can be specified).

## lava 1.5.1

CRAN release: 2017-09-27

- conformal predictions: `confpred`
- warnings (char2num used instead of coersion via as.numeric)
- `%++%` for function compositon
- New `summary.effects` methods with mediation proportion in the output
- New hook: `remove.hooks` (see example `ordinal.lvm`)
- constrain methods now handled more robustly in `sim.lvm` allowing both
  vectorized and non-vectorized functions
- Non-linear associations can now be specified with `nonlinear` method.
  Estimation via the `twostage` function.
- Robust standard errors added to the IV estimator (2SLS)
- New cross-validation function: `cv` (and `csplit` function for
  creating random sets).

## lava 1.5

CRAN release: 2017-03-16

- lava.tobit is longer required for ordinal and censored responses.
  Default is now to use the implementation in the ‘mets’ package.
- Composite likelihood method (`complik`) updated
- weight argument renamed to weights in agreement with lm, glm, coxph, …
- `sim.default`: new argument ‘arg’ passed on to simulation function
- `sim.default`: new argument ‘iter’. If TRUE the iteration number is
  passed to function call as first argument (default FALSE)
- `estimate.default`: Wildcards/global expressions can now be used for
  specifying contrasts based on the syntax of the functions `contr`,
  `parsedesign`. See examples on the help-page. The argument
  transform.ci has been renamed to back.transform.
- correlation methods for matrices and data.frames (either pairwise or
  full MLE). All methods can now return the influence functions.
- `revdiag`: dimnames are kept
- `Combine`: output updated
- `forestplot`: point estimates shown by default
- `backdoor` now works without conditioning set (yields all possible
  conditioning sets)
- New formula syntax: y+x~v+z same as c(y,x)~v+z
- `spaghetti`: trend.formula can now contain a factor statement on the
  rhs

## lava 1.4.7

CRAN release: 2017-01-27

- Maintenance release
- models can now be specified as y1+y2~x1+x2 instead of c(y1,2y)~x1+x2
- `sim` method now has a seed argument

## lava 1.4.6

CRAN release: 2016-12-20

- New backtrace algorithms for Newton-Raphson optimization routine
  (`NR`).
- `diagtest` updated.

## lava 1.4.5

CRAN release: 2016-10-26

- New graph functions: `dsep`: check for d-separation (conditional
  independence). `backdoor`: check backdoor criterion of a graph
  (lvm-object). `adjMat`: return adjaceny matrix. `edgeList`: return
  edge list. `ancestors`: return ancenstors of nodes. `descendants`:
  return descendants of nodes.
- All simple paths in a graph can now be extracted with:
  `path(...,all=TRUE)`
- Covariance parameters are now reference with `~~` instead of `,`.
  Applies to setting starting values in `estimate`, parameters in
  `sim`,`compare`,`estimate`,… To use the old syntax set
  `lava.options(symbol=c("~",","))`.
- `layout` argument added to `lava.options` (default ‘dot’)
- visNetwork support, new `plot.engine` argument added to plot methods.
- `bootstrap.lvmfit` now default returns original estimates.
- `print`, `transform` methods updated (transform output).
- `+` operator overloaded for lvm and estimate objects (merge).
- New composite likelihood function: `complik`.
- New functions for simple association measures: `riskcomp`, `rdiff`,
  `rratio`, …
- New argument ‘latent’ in `simulate` method. If FALSE the latent
  variables are dropped from the returned data.frame.
- `modelsearch` by default now shows both directional or undirectional
  associations (type=‘all’ vs type=‘cor’).
- `sim.default` now stores timings. New print functions (data.table like
  output).
- lvm model can now be updated with the `sim` function, for instance
  setting parameter values for the simulation only once:
  `m <- sim(m,p=p,...)`, with faster subsequent calls `sim(m,n=n)`.
- `estimate.default` can now simulate p-values (‘R’ argument). Returns
  an object which can also be used as input for `estimate`.
- Bug fixes: `NR` optimization with back-tracing; fixed `matrices.lvm`
  when called without variance parameters; fixed a bug in r-square
  computations.
- Contrast matrix can be specified with the function `contr`.

## lava 1.4.4

CRAN release: 2016-08-13

- estimate.default will now use the id-variable of an ‘estimate’ object
  if the ‘id’ argument is left unspecified. For
  multinomial,gkgamma,kappa additional arguments (…) are now parsed on
  the ‘estimate.default’ (including id).
- Updated print/summary methods for ‘estimate.default’.
  Sample/cluster-size added to output.
- Code clean-up and optimization. Smarter calculations of kronecker
  products, and some regular expressions updated.
- New function ‘predictlvm’ which return jacobian.
- Intercepts can now be specified via parantheses, e.g., y ~ (-2) + x
- ‘getoutcome’ with sep argument for splitting ‘\|’ statements in
  formulas.
- Partial gamma, gkgamma, updated (probability interpretation,
  homogeneity tests removed)
- ‘moments’ function now returns conditional mean with multiple rows.
  Side effect fixed across multiple functions
- twostage function with support for mixture models
- Beta (Beta.lvm) and Finite Gaussian (GM2.lvm,GM3.lvm) Mixtures added.
- ‘sim’: parameters can now be specified as part of ‘…’
- summary.sim: calculate Wald CI if confint=TRUE, otherwise use the user
  supplied confidence limits.
- Clopper-pearson intervals and exact binomial tests added to
  ‘diagtest’.
- Interval censoring with ‘normal’ estimator, which now also works with
  ‘binary’ definitions.
- default plot style updated.

## lava 1.4.3

CRAN release: 2016-04-11

- partial gamma coefficients (gkgamma)
- Unit tests works with new testthat version
- Avoid trying to fork new processes on windows (bootstrap,sim.default)

## lava 1.4.2

CRAN release: 2016-04-05

- Code optimization and minor bug fixes
- Travis-CI, unit-tests
- glm estimator update (censored regression)
- polychoric correlations (pcor)
- New utility functions: wrapvec, offdiag
- simulation: regression design on parameters (see weibull + variance
  hetereogeneity example in help(‘sim’))
- Byte compile by default

## lava 1.4.1

CRAN release: 2015-06-22

- New plot.estimate method
- Documentation and examples updated

## lava 1.4.0

CRAN release: 2015-02-17

- Linear measurement error model: ‘measurement.error’
- Diagnostic tests: ‘diagtest’
- ‘plotConf’ updated with support for special function terms (I, poly,
  ns, …). Old version is available (not in namespace) as
  lava:::plotConf0
- Pareto distribution: ‘pareto.lvm’
- Code clean-up/optimization: ‘EventTime’, ‘stack’
- ‘estimate.default’ new syntax for contrast specification (parsedesign)
- ‘regression.lvm’ with y,x argument (as alias for to,from)
- plot longitudinal data: ‘spaghetti’
- Examples updated

## lava 1.3

CRAN release: 2014-11-18

- New syntax for categorical predictors (method ‘categorical’ and
  argument ‘additive=FALSE’ with ’regression method)
- Argument ‘intervals’ added to ‘ones.lvm’ for piece-wise constant
  effects
- Argument ‘average=TRUE’ now needed for empirical averages in
  estimate.default
- Fixed a bug in score.glm (with weights and offset) introduced in
  version 1.2.6
- estimate.default:
  - small-sample corrections
  - Default id from row names in estimate.default (used with merge
    method)
  - iid decompostion also returned for hypothesis contrasts
  - keep argument added to estimate.default and merge
  - labels argument added to estimate.default
- ‘images’ function for visualization of tabular data added to namespace
- ‘ksmooth’ and ‘surface’ for surface estimation and visualization of
  bivariate data and functions
- ‘dsort’: Sort data.frames
- general multivariate distributions in simulations. see example in
  ‘sim’
- ‘or2prob’, ‘tetrachoric’ for conversion from OR to probabilities (and
  tetrachoric correlations). ‘prob.normal’: calculates probabilities
  from threshold model given thresholds and variance See also
  mets:::assoc for calculations of kappa, gamma, uncer.coef.
  ‘normal.threshold’: returns thresholds,variance,mu from model with
  categorical outcomes.
- Multiple testing routines: closed.testing, p.correct, …
- ‘Missing’ method updated with a simple ‘suffix’ argument
- Back-tracing updated in Newton-Raphson routine

## lava 1.2.6

CRAN release: 2014-05-07

- New ‘stack’ function for two-stage estimation (via ‘estimate’ objects)
- New ‘blocksample’ function for resampling clustered data.
- New function ‘Missing’ to generate complex missing data patterns
- Weibull parametrization of ‘coxWeibull.lvm’ rolled back (ver. 1.2.4).
  The function ‘weibull.lvm’ now leads to Accelerated Failure Time model
  (see examples of ‘eventTime’)
- iid function cleanup (new ‘bread’ attribute). iid.glm now gives
  correct estimated influence functions for ‘quasi’ link (constant
  variance)
- Parameter constraints on (co)variance parameters now possible with the
  syntax lvm(…,y\~~a\*x) (corresponding to covariance(…,y~x)\<-“a”)
- Some additional utilities: pdfconvert, scheffe, images, click.
  confband updated with ‘polygon’ argument.
- New function getMplus: Import results from Mplus
- New function getSAS: Import SAS ODS
- New ‘edgecolor’ argument of plot-function

## lava 1.2.5

CRAN release: 2014-03-14

- ‘merge’ method added for combining ‘estimate’ objects
- Adjustments to starting values
- Function ‘categorical’ for adding categorical predictors to simulation
  model
- Improved flexibility in simulations with ‘transform’,‘constrain’ (ex:
  categorical predictors)
- Added ‘dataid’ argument to estimate.default allowing different id for
  ‘data’ and i.i.d. decomposition of model parameter estimates. With the
  argument ‘stack=FALSE’ influence functions within clusters will not be
  stacked together.
- R-squared values (+ approximate standard errors/i.i.d. decomposition)
  via ‘rsq(model,TRUE)’
- New infrastructure for adding additional parameters to models (no
  user-visible changes).
- multinomial function for calculating influence curves for multinomial
  probabilities. ‘gammagk’ and ‘kappa’ methods for calculating
  Goodman-Kruskals gamma and Cohens kappa coefficients.
- ordreg function for univariate ordinal regression models
- iid methods for data.frames/matrices (empirical mean and variance)
- Support for using ‘mets::cluster.index’ in GEE-type models (much
  faster).
- plotConf updated (vcov argument added and more graphical arguments
  parsed to plotting functions)
- Additional unit-tests implemented
- New ‘forestplot’ and ‘Combine’ functions
- Covariance structure may now be specified using ‘~~’, e.g.
  ’lvm(c(y,v)~~z+u)’ specifies correlation between residuals of
  (y,z),(y,u),(v,z),(v,u).

## lava 1.2.4

CRAN release: 2013-12-17

- Avoid estimating IC in ‘estimate.default’ when ‘vcov’ argument is
  given.
- New default starting values
- Time-varying effects via ‘timedep’
- R-squared added to summary
- alias: covariance-\>variance
- added size argument to binomial.lvm;

## lava 1.2.3

CRAN release: 2013-10-28

- ‘subset’ argument added to estimate.default. Calculates empirical
  averages conditional on subsets of data
- Improved output from compare/estimate functions
- Minor bug fixes (plot, predict)
- sim: Piecewise constant rates with coxEponential.lvm. New
  aalenExponential.lvm function for additive models. Functions ones.lvm
  and sequence.lvm for deterministic variables.

## lava 1.2.2

CRAN release: 2013-07-12

- Regression parameters are now by default referenced using ‘~’,
  e.g. “y~x” instead of “y\<-x”. Applies to setting starting values in
  ‘estimate’, parameters in ‘sim’,‘compare’,‘estimate’,…. To use the old
  syntax set ‘lava.options(symbol=c(“\<-”,“\<-\>”))’
- Newton-Raphson/scoring procedure updated
- Search-interval for profile likelihood CI improved (for variance
  parameters)
- ‘estimate.default’ updated (LRT)
- ‘iid’ updated (variance now obtained as tensor product of the result)
- progress bar for ‘bootstrap’ and ‘modelsearch’
- various minor bug fixes
- new functions: Expand (expand.grid wrapper), By (by wrapper)

## lava 1.2.1

CRAN release: 2013-05-21

- Optimization + minor bug fixes

## lava 1.2.0

CRAN release: 2013-04-24

- New method ‘iid’ for extracting i.i.d. decomposition (influence
  functions) from model objects (e.g. glm, lvm, …)
- Method ‘estimate’ can now be used on model objects to transform
  parameters (Delta method) or conduct Wald tests. Average effects,
  i.e. averaging functionals over the empirical distribution is also
  possible including calculation of standard errors.
- ‘curereg’ function for estimating mixtures of binary data.
- Instrumental Variable (IV) estimator (two-stage least-squares)
  optimized.
- New distributions: Gamma.lvm, coxWeibull.lvm, coxExponential.lvm,
  coxGompertz.lvm. New method ‘eventTime’ (for simulation of competing
  risks data)
