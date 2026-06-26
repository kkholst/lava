#' @export
estimate <- function(x, ...) UseMethod("estimate")

#' Influence function based inference
#'
#' Primary tool for obtaining parameter estimates with robust (sandwich)
#' standard errors, applying the delta method, and testing linear hypotheses.
#' The function returns an object of class `estimate` which serves as a
#' general container for parameter estimates and their influence functions
#' (IFs). Three calling conventions are supported:
#'
#' - `estimate(x, ...)` -- extract estimates from a model object
#' - `estimate(coef=, IC=, ...)` -- construct from coefficients and IF matrix
#' - `estimate(coef=, vcov=, ...)` -- construct from coefficients and
#'   covariance matrix
#'
#' @param x model object (`glm`, `lvmfit`, ...) or an existing `estimate`
#'   object. When two model objects are supplied (e.g., `estimate(g, g0)`) a
#'   likelihood-ratio test is performed.
#' @param f transformation of model parameters. Accepts several input types:
#'   - A **function** `f(p)` or `f(p, data)`: applies the delta method.
#'     When `f` returns a named list the names are used as parameter labels.
#'   - A **matrix**: used as a contrast (linear combination) matrix.
#'   - A **numeric vector** of parameter indices: converted to a contrast
#'     that selects and differences those parameters.
#'   - A **list** of indices: each element selects one parameter.
#'   - **Character** expressions: supports wildcards (`"?"`, `"*"`) and
#'     arithmetic on parameter names (e.g., `"z" - "x"`, `2 * "z" - 3 * "x"`).
#' @param ... additional arguments to lower level functions
#' @param data `data.frame` used by `f` when the transformation depends on
#'   covariates (see `average`). Defaults to `model.frame(x)`.
#' @param id (optional) cluster identifier. Can be a vector of cluster IDs, a
#'   one-sided formula (evaluated in `data`), a single character column name,
#'   or a logical scalar (`TRUE` for one-to-one matching, `FALSE` for
#'   independence). When supplied, the IF is aggregated within clusters to
#'   produce cluster-robust standard errors.
#' @param coef (optional) named parameter vector. Used instead of
#'   `coef(x)` when constructing an `estimate` object without a model.
#' @param IC if `TRUE` (default) the influence function matrix is estimated and
#'   stored in the returned object (extract with the [IC] method). Can also be
#'   a user-supplied IF matrix (one row per observation, one column per
#'   parameter), which is used directly instead of estimating it from `x`.
#' @param vcov (optional) covariance matrix of parameter estimates, or a
#'   logical. If `TRUE`, [stats::vcov] is used to obtain the (model-based)
#'   covariance matrix from `x`, yielding non-robust standard errors. If a
#'   matrix is supplied it is used directly. When omitted or `FALSE`, robust
#'   standard errors are computed from the influence function.
#' @param stack if `TRUE` (default) the influence function contributions are
#'   summed within each cluster defined by `id`. Set to `FALSE` to keep the
#'   un-stacked (per-observation) decomposition.
#' @param average if `TRUE` the function computes the standardized
#'   (marginalized) estimate \eqn{\hat\Psi = P_n f(X; \hat\theta)}, i.e.,
#'   the empirical mean of `f(p, data)` over all rows of `data`. The
#'   influence function accounts for both the empirical averaging and the
#'   parameter estimation uncertainty (see Details).
#' @param subset (optional) logical vector, expression evaluated in `data`, or
#'   column name. When used together with `average = TRUE`, the average is
#'   conditioned on the subpopulation where `subset` is `TRUE`, yielding a
#'   conditional marginalized estimate.
#' @param keep (optional) index of parameters to keep from final result.
#'   Accepts integer indices, character names, or (with `regex = TRUE`)
#'   perl-compatible regular expressions.
#' @param use (optional) index of parameters to use in calculations. The
#'   selected parameters are first extracted (via `keep`) and then the
#'   remaining arguments (`f`, `contrast`, etc.) are applied to this subset.
#' @param regex if `TRUE` use perl-compatible regular expressions for `keep`
#'   and `use` arguments
#' @param ignore.case ignore case in regular expressions
#' @param print (optional) custom print function for the resulting `estimate`
#'   object
#' @param labels (optional) character vector of coefficient names
#' @param label.width (optional) max display width of labels
#' @param contrast (optional) contrast matrix for a final Wald test. When
#'   supplied together with `null`, tests \eqn{H_0: B\theta = b_0}.
#' @param null (optional) null hypothesis vector \eqn{b_0} to test against
#'   (default 0)
#' @param level level of confidence limits (default 0.95)
#' @param type type of small-sample correction for cluster-robust variance.
#'   One of:
#'   - `"robust"` (default): no correction.
#'   - `"df"`: applies \eqn{n/(n-p)} correction (Mancl & DeRouen, 2001).
#'   - `"mbn"`: Morel-Bokossa-Neerchal (2003) correction.
#'   - `"hc3"`: leverage-adjusted HC3-type correction (blended with
#'     `var.adj`).
#'   - `"hc4"`: Cribari-Neto (2004) leverage-adjusted correction.
#' @param var.adj blending parameter for the HC3 leverage adjustment
#'   (default 0.25). Controls the weight between observation-level empirical
#'   leverage and the average leverage \eqn{p/n}.
#' @param df degrees of freedom for t-based inference (default: `NULL` for
#'   Gaussian approximation; when set, confidence intervals and p-values use
#'   the t-distribution with `df` degrees of freedom)
#' @param back.transform (optional) function applied to the point estimates
#'   and confidence interval bounds *after* inference is performed on the
#'   original scale. Useful for variance-stabilizing transformations, e.g.,
#'   compute CIs on the `atanh` (Fisher z) scale and back-transform with
#'   `tanh`.
#' @details
#'
#' # Influence functions and robust standard errors
#'
#' An estimator \eqn{\widehat{\theta}} is *regular and asymptotically
#' linear* (RAL) when it admits the iid decomposition
#' \deqn{\sqrt{n}(\widehat{\theta}-\theta) =
#' \frac{1}{\sqrt{n}}\sum_{i=1}^n \mathrm{IC}(Z_i; P) + o_p(1)}
#' where \eqn{\mathrm{IC}} is the unique *influence function* satisfying
#' \eqn{E\{\mathrm{IC}(Z; P)\} = 0}. By the central limit theorem
#' \deqn{\sqrt{n}(\widehat{\theta}-\theta)
#' \overset{d}{\longrightarrow}
#' N(0,\; \mathrm{Var}\{\mathrm{IC}(Z; P)\})}
#' and the asymptotic variance is consistently estimated by the empirical
#' variance of the plugin IF estimate, yielding robust (sandwich) standard
#' errors. The estimated IF can be extracted with the [IC] method.
#'
#' # Parameter transformations (delta method)
#'
#' When `f` is a function \eqn{\phi: R^p \to R^m}, the delta method is
#' applied:
#' \deqn{\sqrt{n}\{\phi(\widehat{\theta}) - \phi(\theta)\} =
#' \frac{1}{\sqrt{n}}\sum_{i=1}^n
#' \nabla\phi(\theta)\,\mathrm{IC}(Z_i; P) + o_p(1)}
#' Derivatives are computed numerically via [numDeriv::jacobian] unless the
#' function returns an attribute `"grad"` with the analytic Jacobian.
#'
#' Alternatively, `estimate` objects support direct arithmetic operations
#' (e.g., `a * b`, `exp(a)`, `a^b`) which apply the delta method with
#' *exact* (analytical) derivatives computed automatically. This influence
#' function calculus allows building complex transformations from simple
#' building blocks without numerical differentiation. See the last example
#' section ("influence function calculus") and
#' `vignette("influencefunction", package = "lava")` for details.
#'
#' # Averaging and marginalization
#'
#' When `average = TRUE` and `f(p, data)` depends on covariates, the
#' target parameter is the standardized (marginalized) estimate
#' \eqn{\Psi = E\{f(X;\theta)\}}. The IF for the averaged estimate
#' accounts for both the empirical averaging and parameter estimation
#' uncertainty:
#' \deqn{\mathrm{IC}_\Psi(Z; P) = f(X;\theta) - \Psi +
#' [E\nabla_\theta f(X;\theta)]\,\phi(Z; P)}
#' When `subset` is also specified, the average is conditioned on the
#' subpopulation, yielding a conditional marginalized estimate.
#'
#' # Cluster-robust standard errors
#'
#' When `id` is supplied, the per-observation IF contributions are summed
#' within clusters (when `stack = TRUE`), producing the cluster-level IF
#' \eqn{\widetilde{\mathrm{IC}}(Z_i; P) = \sum_{k=1}^{N_i}
#' \frac{n}{N}\mathrm{IC}(Z_{ik}; P)}.
#' The resulting variance estimate is equivalent to the GEE working
#' independence sandwich estimator.
#'
#' For full theoretical background and worked examples see
#' `vignette("influencefunction", package = "lava")`.
#'
#' @export
#' @export estimate.default
#' @examples
#'
#' ## Simulation from logistic regression model
#' m <- lvm(y~x+z);
#' distribution(m,y~x) <- binomial.lvm("logit")
#' d <- sim(m,1000)
#' g <- glm(y~z+x,data=d,family=binomial())
#' g0 <- glm(y~1,data=d,family=binomial())
#'
#' ## LRT
#' estimate(g, g0)
#'
#'
## Plain estimates (robust standard errors)
#' estimate(g)
#'
#' ## Testing contrasts
#' summary(estimate(g), null=0)
#' estimate(g, rbind(c(1,1,0), c(1,0,2)))
#' summary(estimate(g, rbind(c(1,1,0), c(1,0,2))), null=c(1,2))
#' estimate(g, 2:3) ## same as cbind(0,1,-1)
#' estimate(g, as.list(2:3)) ## same as rbind(c(0,1,0),c(0,0,1))
#' ## Alternative syntax
#' estimate(g, "z", "z"-"x", 2*"z"-3*"x")
#' estimate(g, "?")  ## Wildcards
#' estimate(g, "*Int*", "z")
#' estimate(g, "1", "2"-"3", null = c(0,1))
#' estimate(g, 2, 3)
#'
#' ## Usual (non-robust) confidence intervals
#' estimate(g, vcov=TRUE)
#' estimate(g, vcov=vcov(g))
#'
#' ## Transformations
#' estimate(g, function(p) p[1]+p[2])
#'
#' ## Multiple parameters
#' e <- estimate(g, function(p) c(p[1]+p[2], p[1]*p[2]))
#' e
#' vcov(e)
#'
#' ## Label new parameters
#' estimate(g, function(p) list("a1"=p[1]+p[2], "b1"=p[1]*p[2]))
#' #'
#' ## Multiple group
#' m <- lvm(y~x)
#' m <- baptize(m)
#' d2 <- d1 <- sim(m,50,seed=1)
#' e <- estimate(list(m,m),list(d1,d2))
#' estimate(e) ## Wrong
#' ee <- estimate(e, id=rep(seq(nrow(d1)), 2)) ## Clustered
#' ee
#' estimate(lm(y~x,d1))
#'
#' ## Marginalize
#' f <- function(p,data)
#'   list(p0=expit(p["(Intercept)"] + p["z"]*data[,"z"]),
#'        p1=expit(p["(Intercept)"] + p["x"] + p["z"]*data[,"z"]))
#' e <- estimate(g, f, average=TRUE)
#' e
#' estimate(e,diff)
#' estimate(e,cbind(1,1))
#'
#' ## Clusters and subset (conditional marginal effects)
#' d$id <- rep(seq(nrow(d)/4),each=4)
#' estimate(g,function(p,data)
#'          list(p0=expit(p[1] + p["z"]*data[,"z"])),
#'          subset=d$z>0, id=d$id, average=TRUE)
#'
#' ## More examples with clusters:
#' m <- lvm(c(y1,y2,y3)~u+x)
#' d <- sim(m,10)
#' l1 <- glm(y1~x,data=d)
#' l2 <- glm(y2~x,data=d)
#' l3 <- glm(y3~x,data=d)
#'
#' ## Some random id-numbers
#' id1 <- c(1,1,4,1,3,1,2,3,4,5)
#' id2 <- c(1,2,3,4,5,6,7,8,1,1)
#' id3 <- seq(10)
#'
#' ## Un-stacked and stacked i.i.d. decomposition
#' IC(estimate(l1,id=id1,stack=FALSE))
#' IC(estimate(l1,id=id1))
#'
#' ## Combined i.i.d. decomposition
#' e1 <- estimate(l1,id=id1)
#' e2 <- estimate(l2,id=id2)
#' e3 <- estimate(l3,id=id3)
#' (a2 <- merge(e1,e2,e3))
#'
#' ## If all models were estimated on the same data we could use the
#' ## syntax:
#' ## Reduce(merge,estimate(list(l1,l2,l3)))
#'
#' ## Same:
#' IC(a1 <- merge(l1,l2,l3,id=list(id1,id2,id3)))
#'
#' IC(merge(l1,l2,l3,id=TRUE)) # one-to-one (same clusters)
#' IC(merge(l1,l2,l3,id=FALSE)) # independence
#'
#'
#' # ------ influence function calculus -------
#' a <- estimate(coef = c("a" = 0.5), IC = scale(rnorm(10), scale=FALSE), id = 1:10)
#' b <- estimate(coef = c("b" = 0.8), IC = scale(rnorm(10), scale=FALSE), id = 1:10)
#'
#' e <- c(a, b) # merge
#' merge(a, b)
#' c(e1=a, b) # naming of par
#' labels(e, c("p1", "p2")) # renaming parameters
#' e["a"] # subset
#' subset(e, "a")
#'
#' # pipes
#' # c(a, b) |>
#' #  transform(function(x) x^2) |>
#' #  subset("a") |>
#' #  labels("sq")
#'
#' # Parameter transformation with automatic calculation of derivatives
#' a * b
#' (3 * cos(a) / sqrt(b) + 1) / a
#' expit(c(a,b))
#' c(sum=sum(e), sum2=a+b,
#'   prod=prod(e), prod2=a*b)
#' e %*% e # inner prod.
#' c(1, 2) %*% e
#' c(pow = a^b)
#' a^c(0.5, 2)
#' c(b=e["a"] * e["b"] / a, also.b=e["b"])
#'
#' B <- rbind(c(1,-1), c(1,0), c(0,1))
#' B %*% e
#' e == 1 # wald-test, null-hypothesis H0: b=1
#' e == c(1,2)
#' B %*% e == 1
#' @aliases estimate estimate.default
#' @aliases estimate.mlm
#' @seealso [estimate.array], [merge.estimate], [contr], [parsedesign],
#'   [pairwise.diff], [c.estimate], [summary.estimate],
#'   `coef.estimate`,
#'   `vcov.estimate`, `transform.estimate`, `labels.estimate`,
#' @return Object of class `estimate` with the following elements:
#'   \item{coef}{Named vector of parameter estimates.}
#'   \item{vcov}{Variance-covariance matrix.}
#'   \item{IC}{Influence function matrix (observations x parameters).}
#'   \item{coefmat}{Formatted coefficient table (estimate, std.err,
#'     confidence limits, p-value).}
#'   \item{id}{Cluster/id variable used.}
#'   \item{ncluster}{Number of clusters.}
#'   \item{n}{Number of observations.}
#'   \item{compare}{(When `null` or contrasts are specified) Wald test result.}
#' @method estimate default
#' @export
estimate.default <- function(x=NULL, f=NULL, ...,
                             data, id,
                             coef, IC=TRUE, vcov,
                             stack=TRUE,
                             average=FALSE, subset,
                             keep, use,
                             regex=FALSE, ignore.case=FALSE,
                             print=NULL, labels, label.width,
                             # deprecated argumets (removed in version 1.9.3)
                             # return object of class summary.estimate
                             contrast,
                             null,
                             level=NULL,
                             type=NULL,
                             var.adj=NULL,
                             df=NULL,
                             back.transform=NULL
                             ) {
  cl <- match.call(expand.dots = TRUE)
  cal <- match.call()

  if ("iid" %in% names(cl)) {
    stop("The 'iid' argument is obsolete. Please use the 'IC' argument")
  }
  if ("R" %in% names(cl) || "null.sim" %in% names(cl)) {
    stop("The 'R' and 'null.sim' arguments have been removed. ")
  }
  if ("robust" %in% names(cl)) {
    warning(
      "The 'robust' argument is deprecated and ignored. ",
      "Robust standard errors are now always computed. ",
      "Use the 'vcov' argument to compute model-based SEs."
    )
  }

  if (!missing(use)) {
    p0 <- c(
      "f", "subset", "average",
      "keep", "labels", "null"
    )
    cl0 <- cl
    cl0[c("use", p0)] <- NULL
    cl0$keep <- use
    cl$x <- eval(cl0, parent.frame())
    cl[c("vcov", "use")] <- NULL
    return(eval(cl, parent.frame()))
  }
  expr <- suppressWarnings(inherits(try(f, silent=TRUE), "try-error"))
  if (!missing(coef)) {
    pp <- coef
  } else {
    pp <- suppressWarnings(try(stats::coef(x), silent = TRUE))
    if (inherits(x, "survreg") && length(pp) < NROW(x$var)) {
      pp <- c(pp, scale=x$scale)
    }
  }
  if (expr || is.character(f) || (is.numeric(f)
    && !is.matrix(f))) { ## || is.call(f)) {
    dots <- lapply(substitute(placeholder(...))[-1], function(x) x)
    args <- c(list(
      coef = names(pp),
      x = substitute(f),
      regex = regex
    ), dots)
    f <- do.call(parsedesign, args)
  }

  contrast.transform <- FALSE  # if TRUE parameter estimates should be
                               # transformed according to contrast matrix 'f'
  if (!is.null(f) && !is.function(f)) {
    if (!(is.matrix(f) || is.vector(f)))
      return(compare(x, f, ...)) ## LRT
    contrast.transform <- TRUE
  }

  if (lava.options()$cluster.index) {
    if (!requireNamespace("mets", quietly=TRUE)) stop("'mets' package required")
  }
  if (missing(data))
    data <- tryCatch(model.frame(x), error=function(...) NULL)
  nn <- NULL
  if (((is.logical(IC) && IC) || length(IC)>0) &&
      (missing(vcov) || is.null(vcov) ||
       (is.logical(vcov) && vcov[1]==FALSE && !is.na(vcov[1])))) {
    ## If user supplied vcov, then don't estimate IC
    if (!is.logical(IC)) {
      ic_theta <- cbind(IC)
      if (NCOL(ic_theta) != length(pp)) {
        warning("Wrong dimension of influence function IC")
      }
      if (lava.options()$check.ic) {
        check_ic_mean_zero(ic_theta)
      }
      IC <- TRUE
    } else {
      ic_theta <- IC(x)
    }
  } else {
    if (!is.null(x) && (missing(vcov) ||
                        (is.logical(vcov) && !is.na(vcov)[1])))
      suppressWarnings(vcov <- stats::vcov(x))
    ic_theta <- NULL
  }

  if (any(is.na(ic_theta))) {
    ## Rescale each column according to I(obs)/pr(obs)
    for (i in seq_len(NCOL(ic_theta))) {
      pr <- mean(!is.na(ic_theta[, i]))
      ic_theta[, i] <- ic_theta[, i]/pr
    }
    ic_theta[is.na(ic_theta)] <- 0
  }

  if (!missing(subset)) {
    e <- substitute(subset)
    expr <- suppressWarnings(inherits(try(subset, silent=TRUE), "try-error"))
    if (expr) subset <- eval(e, envir=data)
    if (is.character(subset)) subset <- data[, subset]
    if (is.numeric(subset)) subset <- subset > 0
  }
  idstack <- NULL
  ## Preserve id from 'estimate' object
  if (missing(id)) {
    if (inherits(x, "measurement.error")) {
      if (!is.null(x[["id"]])) id <- x[["id"]]
    } else if (inherits(x, "estimate") && !is.null(index(x))) {
      id <- index(x)
    }
  }
  if (!missing(id) && IC) {
    if (is.null(ic_theta)) stop("'IC' method needed")
    nprev <- nrow(ic_theta)
    if (inherits(id, "formula")) {
      id <- interaction(get_all_vars(id, data))
    }
    ## e <- substitute(id)
    ## expr <- suppressWarnings(inherits(try(id,silent=TRUE),"try-error"))
    ## if (expr) id <- eval(e,envir=data)
    ##if (!is.null(data)) id <- eval(e, data)
    if (is.logical(id) && length(id)==1) {
      id <- if(is.null(ic_theta)) seq_len(nrow(data)) else seq_len(nprev)
      stack <- FALSE
    }
    if (is.character(id) && length(id)==1)
      id <- data[, id, drop=TRUE]
    if (!is.null(ic_theta)) {
      if (length(id)!=nprev) {
        if (!is.null(x$na.action) &&
            (length(id)==length(x$na.action) + nprev)) {
          warning("Applying na.action")
          id <- id[-x$na.action]
        } else stop("Dimensions of i.i.d decomposition and 'id' does not agree")
      }
    } else {
      if (length(id)!=nrow(data)) {
        if (!is.null(x$na.action) &&
            (length(id)==length(x$na.action)+nrow(data))) {
          warning("Applying na.action")
          id <- id[-x$na.action]
        } else stop("Dimensions of IC and 'id' does not agree")
      }
    }
    if (stack) {
      N <- nrow(ic_theta)
      clidx <- NULL
      atr <- attributes(ic_theta)
      atr$dimnames <- NULL
      atr$dim <- NULL
      if (!lava.options()$cluster.index) {
        ic_theta <- matrix(unlist(by(ic_theta, id, colSums)),
                           byrow=TRUE, ncol=ncol(ic_theta))
        attributes(ic_theta)[names(atr)] <- atr
        idstack <- sort(unique(id))
      } else {
        clidx <- mets::cluster.index(id, mat=ic_theta, return.all=TRUE)
        ic_theta <- with(clidx, X)
        attributes(ic_theta)[names(atr)] <- atr
        idstack <- id[as.vector(clidx$firstclustid)+1]
      }
      ic_theta <- ic_theta*NROW(ic_theta)/length(id)
      if (is.null(attributes(ic_theta)$N)) {
        attributes(ic_theta)$N <- N
      }
    } else idstack <- id
  } else {
    if (!is.null(data)) idstack <- rownames(data)
  }
  if (!is.null(ic_theta) && (length(idstack)==nrow(ic_theta))) {
    rownames(ic_theta) <- idstack
  }
  if (!is.null(ic_theta) && (missing(vcov) || is.null(vcov))) {
    V <- var_ic(ic_theta)
  } else {
    if (!missing(vcov)) {
      if (length(vcov) == 1 && is.na(vcov)) {
        vcov <- matrix(NA, length(pp), length(pp))
      }
      V <- cbind(vcov)
    } else {
      suppressWarnings(V <- stats::vcov(x))
    }
  }

  if (contrast.transform) {
    B <- f
    if (is.vector(B) || is.list(B)) {
      B <- contr(f, names(pp), ...)
    }
    obj <- structure(list(coef=pp, vcov=V), class="estimate")
    cc <- compare(obj, contrast=B) # to construct new parameter names
    pp <- as.vector(B %*% pp)
    names(pp) <- strip_bracket(cc$cnames)
    if (!is.null(ic_theta)) {
      ic_theta <- ic_theta %*% t(B)
    }
    V <- B %*% V %*% t(B)
    f <- NULL
  }

  derivative <- NULL
  if (!is.null(f)) {
    form <- names(formals(f))
    dots <- ("..."%in%names(form))
    form0 <- setdiff(form, "...")
    parname <- "p"
    if (!is.null(form)) parname <- form[1] # unless .Primitive
    if (length(form0)==1 && !(form0%in%c("object", "data"))) {
      parname <- form0
    }
    if (!is.null(ic_theta)) {
      arglist <- c(list(object=x, data=data, p=vec(pp)), list(...))
      names(arglist)[3] <- parname
    } else {
      arglist <- c(list(object=x, p=vec(pp)), list(...))
      names(arglist)[2] <- parname
    }
    if (!dots) {
      arglist <- arglist[intersect(form0, names(arglist))]
    }
    newf <- NULL
    if (length(form)==0) {
      arglist <- list(vec(pp))
      newf <- function(...) do.call("f", list(...))
      val <- do.call("f", arglist)
    } else {
      val <- do.call("f", arglist)
      if (is.list(val)) {
        nn <- names(val)
        val <- do.call("cbind", val)
        newf <- function(...) do.call("cbind", f(...))
      }
    }
    k <- NCOL(val)
    N <- NROW(val)
    D <- attributes(val)$grad
    if (!is.null(D)) derivative <- D
    if (is.null(D)) {
      D <- numDeriv::jacobian(function(p, ...) {
        if (length(form)==0) arglist[[1]] <- p
        else arglist[[parname]] <- p
        if (is.null(newf))
          return(do.call("f", arglist))
        return(do.call("newf", arglist)) }, pp)
    }
    if (is.null(ic_theta)) {
      pp <- structure(as.vector(val), names=names(val))
      V <- D%*%V%*%t(D)
    } else {
      if (!average || (N<NROW(data))) {  ## transformation not depending on data
        pp <- structure(as.vector(val), names=names(val))
        ic_theta <- ic_theta%*%t(D)
        V <- var_ic(ic_theta)
      } else {
        if (k>1) { ## More than one parameter (and depends on data)
          if (!missing(subset)) { ## Conditional estimate
            val <- apply(val, 2, function(x) x*subset)
          }
          D0 <- matrix(nrow=k, ncol=length(pp))
          for (i in seq_len(k)) {
            D1 <- D[seq(N)+(i-1)*N, , drop=FALSE]
            if (!missing(subset)) ## Conditional estimate
              D1 <- apply(D1, 2, function(x) x*subset)
            D0[i, ] <- colMeans(D1)
          }
          D <- D0
          ic2 <- ic_theta%*%t(D)
        } else { ## Single parameter
          if (!missing(subset)) { ## Conditional estimate
            val <- val*subset
            D <- apply(rbind(D), 2, function(x) x*subset)
          }
          D <- colMeans(rbind(D))
          ic2 <- ic_theta%*%D
        }
        pp <- vec(colMeans(cbind(val)))
        ic1 <- (cbind(val)-rbind(pp)%x%cbind(rep(1, N)))
        if (NROW(ic_theta)==NROW(ic1)) {
          rownames(ic1) <- rownames(ic_theta)
        }
        if (!missing(id)) {
          if (!lava.options()$cluster.index)
            ic1 <- matrix(unlist(by(ic1, id, colSums)),
                          byrow=TRUE, ncol=ncol(ic1))
          else {
            ic1 <- mets::cluster.index(id, mat=ic1, return.all=FALSE)
          }
          ic1 <- ic1 * NROW(ic1) / length(id)
        }
        if (!missing(subset)) { ## Conditional estimate
          phat <- mean(subset)
          ic3 <- cbind(-1/phat^2 * (subset-phat)) ## check
          if (!missing(id)) {
            if (!lava.options()$cluster.index) {
              ic3 <- matrix(unlist(by(ic3, id, colSums)),
                            byrow=TRUE, ncol=ncol(ic3))
            } else {
              ic3 <- mets::cluster.index(id, mat=ic3, return.all=FALSE)
            }
          }
          ic3 <- ic3*NROW(ic3)/length(id)
          ic_theta <- (ic1+ic2)/phat + rbind(pp)%x%ic3
          pp <- pp/phat
          V <- var_ic(ic_theta)
        } else {
          if (nrow(ic1)!=nrow(ic2)) {
            message("Assuming independence between model iid decomposition and new data frame") #nolint
            V <- var_ic(ic1) + var_ic(ic2)
          } else {
            ic_theta <- ic1+ic2
            V <- var_ic(ic_theta)
          }
        }
      }
    }
  }

  df <- NULL
  if (inherits(x, "lm") && family(x)$family == "gaussian"
      && !missing(vcov)) {
    # defaults to t-distribution when calculating p-values with model-based SEs
    df <- x$df.residual
  }
  if (is.null(V)) {
    res <- cbind(pp, NA, NA, NA, NA)
  } else {
    if (length(pp)==1)
      res <- rbind(c(pp, diag(V)^0.5))
    else
      res <- cbind(pp, diag(V)^0.5)
  }
  res <- estimate_coefmat(res[, 1], res[, 2], df=df, level=0.95, null=0)
  if (nrow(res)>0)
    if (!is.null(nn)) {
      rownames(res) <- nn
    } else {
      nn <- attributes(res)$varnames
      if (!is.null(nn))
        rownames(res) <- nn
      if (is.null(rownames(res)))
        rownames(res) <- paste0("p", seq_len(nrow(res)))
    }

  if (NROW(res)==0L) {
    coefs <- NULL
  } else {
    coefs <- res[, 1, drop=TRUE]
    names(coefs) <- rownames(res)
  }
  res <- structure(list(coef=coefs, coefmat=res, vcov=V,
                        IC=NULL, print=print, id=idstack, df=df),
                   class="estimate")
  if (IC) {
    res$IC <- ic_theta
  }
  if (length(coefs)==0L) return(res)

  if (!missing(keep) && !is.null(keep)) {
    if (is.character(keep)) {
      if (regex) {
        nn <- rownames(res$coefmat)
        keep <- unlist(lapply(keep, function(x) {
          grep(x, nn,
               perl = TRUE,
               ignore.case = ignore.case
               )
        }))
      } else {
        keep <- match(keep, rownames(res$coefmat))
      }
    }
    res$coef <- res$coef[keep]
    res$coefmat <- res$coefmat[keep, , drop=FALSE]
    if (!is.null(res$IC)) res$IC <- res$IC[, keep, drop=FALSE]
    res$vcov <- res$vcov[keep, keep, drop=FALSE]
  }

  res <- labels.estimate(
    object = res, str = labels, label.width = label.width
  )
  res$call <- cal
  res$df <- df
  res$n <- nrow(data)
  res$ncluster <- if (!is.null(ic_theta)) nrow(ic_theta) else nrow(data)
  res$derivative <- derivative
  res <- structure(res, class="estimate")

  if (!missing(null) || !missing(contrast) ||
      !is.null(type) || !is.null(var.adj) ||
      !is.null(back.transform) || !is.null(level) ||
      !is.null(df)
      ) {
    .Deprecated(
      msg = paste0(
        "The 'null', 'contrast', 'type', 'back.transform', 'level'
        and 'var.adj' arguments of ",
        "estimate.default() are deprecated. Use ",
        "summary(estimate(...),
null=, contrast=, type=, transform=, level=, df=, var.adj=) instead."
      )
    )
    args <- list(object=res)
    if (!missing(contrast)) args$contrast <- contrast
    if (!missing(null)) args$null <- null
    if (!is.null(level)) args$level <- level
    if (!is.null(type)) args$type <- type
    if (!is.null(var.adj)) args$var.adj <- var.adj
    if (!is.null(df)) args$df <- df
    if (!is.null(back.transform)) args$transform <- back.transform
    return(do.call(summary, args))
  }

  return(res)
}

#' @export
print.estimate <- function(x, type=0L, digits=4L, width=25L,
                           std.error=TRUE, p.value=TRUE,
                           sep=cli::symbol[["line"]],
                           sep.which,
                           sep.labels=NULL,
                           indent=" ", unique.names=TRUE,
                           na.print="", ...) {

  if (!is.null(x$print)) {
    x$print(x, digits=digits, width=width, ...)
    return(invisible(x))
  }
  if (type>0 && !is.null(x$call)) {
    cat("Call: ")
    print(x$call)
    print(cli::rule(width=min(cli::console_width(),60)))
  }
  if (type>0) {
    if (!is.null(x[["n"]]) && !is.null(x[["k"]])) {
      cat("n = ", x[["n"]], ", clusters = ", x[["k"]], "\n\n", sep="")
    } else {
      if (!is.null(x[["n"]])) {
        cat("n = ", x[["n"]], "\n\n", sep="")
      }
      if (!is.null(x[["k"]])) {
        cat("n = ", x[["k"]], "\n\n", sep="")
      }
    }
  }

  cc <- x$coefmat
  if (!is.null(rownames(cc)) && unique.names)
    rownames(cc) <- make.unique(
      unlist(lapply(rownames(cc),
                    function(x) toString(x, width=width)))
    )
  if (!std.error) cc <- cc[, -2, drop=FALSE]
  if (!p.value) cc <- cc[, -ncol(cc), drop=FALSE]

  sep.pos <- c()
  if (missing(sep.which) && !is.null(x$model.index)) {
    sep.which <- unlist(lapply(x$model.index,
                               function(x)
                                 tail(x, 1)))[-length(x$model.index)]
  }
  if (missing(sep.which)) sep.which <- NULL

  if (!is.null(sep.which)) {
    sep0 <- 0%in%sep.which
    if (sep0)
      sep.which <- setdiff(sep.which, 0)
    cc0 <- c()
    sep.which <- c(0, sep.which, nrow(cc))
    N <- length(sep.which)-1
    for (i in seq(N)) {
      if ((sep.which[i]+1)<=nrow(cc))
        cc0 <- rbind(cc0, cc[seq(sep.which[i]+1, sep.which[i+1]), , drop=FALSE])
      if (i<N) {
        cc0 <- rbind(cc0, NA)
        sep.pos <- c(sep.pos, nrow(cc0))
      }
    }
    if (sep0) {
      sep.pos <- c(1, sep.pos+1)
      cc0 <- rbind(NA, cc0)
    }
    cc <- cc0
  }
  if (!is.null(sep.labels)) {
    sep.labels <- rep(sep.labels, length.out=length(sep.pos))
    rownames(cc)[sep.pos] <- sep.labels
    rownames(cc)[-sep.pos] <- paste0(indent, rownames(cc)[-sep.pos])
  } else {
    if (length(sep.pos)>0)
      rownames(cc)[sep.pos] <- rep(paste0(rep(sep, max(nchar(rownames(cc)))),
                                          collapse=""), length(sep.pos))
  }
  print(cc, digits=digits, na.print=na.print, ...)

  if (!is.null(x$compare)) {
    print(cli::rule(width=min(cli::console_width(),60)))
    cat(x$compare$method[3], "\n")
    cat(paste(" ", x$compare$method[-(1:3)], collapse="\n"), "\n")
    if (length(x$compare$method)>=4) {
      out <- character()
      out <- with(x$compare, c(out, paste(names(statistic),
                                          "=", format(round(statistic, 4)))))
      out <- with(x$compare, c(out, paste(names(parameter),
                                          "=", format(round(parameter, 3)))))
      fp  <- with(x$compare, format.pval(p.value, digits = digits))
      out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "<")
                                       fp else paste("=", fp)))
      cat(" ", strwrap(paste(out, collapse = ", ")), sep = "\n")
    }
  }
}

#' @export
vcov.estimate <- function(object, list=FALSE, ...) {
  res <- object$vcov
  nn <- names(coef(object, ...))
  if (list && !is.null(object$model.index)) {
    return(lapply(object$model.index, function(x) object$vcov[x, x]))
  }
  dimnames(res) <- list(nn, nn)
  res
}

#' @export
coef.estimate <- function(object,
                          mat=FALSE,
                          list=FALSE,
                          ...) {
  if (mat) return(object$coefmat)
  if (list && !is.null(object$model.index)) {
    return(lapply(object$model.index, function(x) object$coef[x]))
  }
  object$coef
}

#' @export
transform.estimate <- function(`_data`, ...) {
  estimate(`_data`, ...)
}

#' @export
labels.estimate <- function(object, str, label.width, ...) {
  if (!missing(str)) {
    names(object$coef) <- str
    if (!is.null(object$IC))
      colnames(object$IC) <- str
    if (!is.null(object$vcov))
      colnames(object$vcov) <- rownames(object$vcov) <- str
    rownames(object$coefmat) <- str
  }
  if (!missing(label.width)) {
    rownames(object$coefmat) <- make.unique(
      unlist(lapply(rownames(object$coefmat),
                    function(x) toString(x, width = label.width)))
    )
  }
  return(object)
}

#' @export
parameter.estimate <- function(x, ...) {
  return(x$coefmat)
}

#' @export
IC.estimate <- function(x, ...) {
  if (is.null(x$IC)) return(NULL)
  dimn <- dimnames(x$IC)
  if (!is.null(dimn)) {
    dimn[[2]] <- names(coef(x))
  } else {
    dimn <- list(NULL, names(coef(x)))
  }
  structure(x$IC, dimnames=dimn)
}

#' @export
model.frame.estimate <- function(formula, ...) {
  NULL
}
