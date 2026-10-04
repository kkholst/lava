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
#'
#' - A **function** `f(p)` or `f(p, data)`: applies the delta method. When `f`
#'   returns a named list the names are used as parameter labels.
#' - A **matrix**: used as a contrast (linear combination) matrix. - A **numeric
#'   vector** of parameter indices: converted to a contrast that selects and
#'   differences those parameters.
#' - A **list** of indices: each element selects one parameter.
#' - **Character** expressions: supports wildcards (`"?"`, `"*"`) and arithmetic
#'   on parameter names (e.g., `"z" - "x"`, `2 * "z" - 3 * "x"`).
#' @param ... additional arguments to lower level functions
#' @param data `data.frame` used by `f` when the transformation depends on
#'   covariates (see `average`). Defaults to `model.frame(x)`.
#' @param id (optional) cluster identifier. Can be a vector of cluster IDs, a
#'   one-sided formula (evaluated in `data`), a single character column name, or
#'   a logical scalar (`TRUE` for one-to-one matching, `FALSE` for
#'   independence). When supplied, the IF is aggregated within clusters to
#'   produce cluster-robust standard errors. When `average = TRUE`, `id`
#'   refers to the rows of `data` (default: `rownames(data)`), and the
#'   ids of an `estimate` object (`index(x)`) are used as is (non-overlapping
#'   ids are treated as independent observations), whereas the rows of other
#'   model objects are linked via the row names of `data`.
#'   `id = NULL` removes the id (index) and the row names of the influence
#'   function from the returned object.
#' @param coef (optional) named parameter vector. Used instead of `coef(x)` when
#'   constructing an `estimate` object without a model.
#' @param IC if `TRUE` (default) the influence function matrix is estimated and
#'   stored in the returned object (extract with the [IC] method). Can also be a
#'   user-supplied IF matrix (one row per observation, one column per
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
#'   (marginalized) estimate \eqn{\hat\Psi = P_n f(X; \hat\theta)}, i.e., the
#'   empirical mean of `f(p, data)`, as defined by the \code{f} argument,
#'   over all rows of `data`. The influence function accounts for both the
#'   empirical averaging and the parameter estimation uncertainty (see Details).
#' @param subset (optional) logical vector, expression evaluated in `data`, or
#'   column name. When used together with `average = TRUE`, the average is
#'   conditioned on the subpopulation where `subset` is `TRUE`, yielding a
#'   conditional marginalized estimate.
#' @param keep (optional) index of parameters to keep from final result. Accepts
#'   integer indices, character names, or (with `regex = TRUE`) perl-compatible
#'   regular expressions.
#' @param use (optional) index of parameters to use in calculations. The
#'   selected parameters are first extracted (via `keep`) and then the remaining
#'   arguments (`f`, `contrast`, etc.) are applied to this subset.
#' @param regex if `TRUE` use perl-compatible regular expressions for `keep` and
#'   `use` arguments
#' @param ignore.case ignore case in regular expressions
#' @param print (optional) custom print function for the resulting `estimate`
#'   object
#' @param labels (optional) character vector of coefficient names
#' @param label.width (optional) max display width of labels
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
#' The model may be estimated on a different (e.g., smaller) dataset than
#' `data`. The two terms of the IF are then aligned by id: each term is zero
#' for ids outside its own support and rescaled by the inverse proportion of
#' observed ids (as in [merge.estimate]). If there are no common ids the
#' model estimate and `data` are treated as independent.
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
#' distribution(m,y~x) <- dist_bernoulli("logit")
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
#' summary(estimate(g, "1", "2"-"3"), null = c(0,1))
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
#' ## Standardization (g-computation)
#' f <- function(p,data)
#'   list(p0=expit(p["(Intercept)"] + p["z"]*data[,"z"]),
#'        p1=expit(p["(Intercept)"] + p["x"] + p["z"]*data[,"z"]))
#' e <- estimate(g, f, average=TRUE)
#' e
#' estimate(e,diff)
#' estimate(e,cbind(1,1))
#'
#' # g-computation on non-overlapping data:
#' d$id <- paste0("i", 1:nrow(d))
#' d$w <- rbinom(nrow(d), 1, 0.5)
#' d1 <- subset(d, w == 1)
#' g1 <- glm(y ~ x + z, data=d1, family=binomial)
#' e1 <- estimate(g1, id=d1$id)
#' estimate(g1, f, data=d, id="id", average=TRUE)
#'
#' ## Clusters and subset (conditional marginal effects)
#' d$id <- rep(seq(nrow(d)/4),each=4)
#' estimate(g,function(p,data)
#'          list(p0=expit(p[1] + p["z"]*data[,"z"])),
#'          subset=d$z>0, id=d$id, average=TRUE)
#'
#' ## Model estimated on a subset, standardized over the full data
#' g1 <- glm(y~z+x, data=subset(d, id<=100), family=binomial())
#' estimate(g1, function(p,data) expit(p[1] + p["z"]*data[,"z"]),
#'          data=d, id=d$id, average=TRUE)
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
#' ic1 <- scale(rnorm(10), scale=FALSE)
#' a <- estimate(coef = c("a" = 0.5), IC = ic1, id = 1:10)
#' b <- estimate(coef = c("b" = 0.8), IC = ic1, id = 1:10)
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
#'   [pairwise_diff], [c.estimate], [summary.estimate],
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
                             print=NULL, labels, label.width
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
    stop(
      "The 'robust' argument is deprecated. ",
      "Robust standard errors are now always computed. ",
      "Use the 'vcov' argument to compute model-based SEs."
    )
  }
  if (any(c("null", "contrast", "type", "var.adj",
            "back.transform", "level", "df") %in% names(cl))) {
    stop(
        "The 'null', 'contrast', 'type', 'back.transform', 'level' ",
        "and 'var.adj' arguments of estimate.default() are deprecated. ",
        "Use ",
        "summary(estimate(...), null=, contrast=, type=, transform=,",
        "level=, df=, var.adj=) instead."
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
    res <- eval(cl, parent.frame())
    res$call <- cal
    return(res)
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

  if (is.null(names(pp))) {
    names(pp) <- paste0("p", seq_along(pp))
  }

  if ((expr || is.character(f) || (is.numeric(f))
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

  if (missing(data))
    data <- tryCatch(model.frame(x), error=function(...) NULL)
  nn <- NULL
  if (
    (
      (is.logical(IC) && IC) || (length(IC)>0 && !is.logical(IC))) &&
    (missing(vcov) ||
      is.null(vcov) ||
      (is.logical(vcov) && vcov[1]==FALSE && !is.na(vcov[1])))
    ) {
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
  id_user <- !missing(id)
  id_drop <- id_user && is.null(id) # id=NULL: remove id (index) from result
  ## Standardization (average=TRUE): 'id' refers to the rows of 'data', and the
  ## model IF is aligned to these ids (see 'average_align_ids')
  avg_align <- isTRUE(average) && is.function(f) &&
    !is.null(ic_theta) && !is.null(data)
  id_data <- NULL
  if (avg_align) {
    ids <- average_align_ids(x, data=data, ic=ic_theta,
                             id=if (id_user) id else NULL)
    id_data <- ids$id_data
    ic_theta <- ids$ic
    idstack <- ids$id_model
  }
  if (!avg_align) {
    ## Cluster id of the rows of the IF (default: id of 'estimate' object)
    id0 <- if (id_user) id
           else if (inherits(x, "measurement.error")) {
               if (!is.null(x[["id"]])) x[["id"]]
             } else if (inherits(x, "estimate")) index(x)
    if (!is.null(id0) && IC) {
      if (is.null(ic_theta)) stop("'IC' method needed")
      n <- nrow(ic_theta)
      if (is.logical(id0) && length(id0) == 1) stack <- FALSE
      id0 <- resolve_id(id0, data, n=n, x=x, default=seq_len(n))
      if (stack) {
        ic_theta <- cluster_sum_ic(ic_theta, id0)
        idstack <- unique(id0)
      } else {
        idstack <- id0
      }
    } else if (!id_drop && !is.null(data)) {
      idstack <- rownames(data)
    }
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
        ## Empirical averaging term, f(X; theta) - Psi, on the ids of 'data'
        ic1 <- (cbind(val)-rbind(pp)%x%cbind(rep(1, N)))
        ic1 <- cluster_sum_ic(ic1, id_data)
        uid_data <- unique(id_data)
        uid_model <- idstack
        if (!any(uid_model %in% uid_data)) {
          message("Assuming independence between model iid decomposition and new data frame") #nolint
        }
        ## Align the terms by id (union of ids, starting with the ids of
        ## 'data'), see 'align_ic'
        if (!missing(subset)) { ## Conditional estimate
          phat <- mean(subset)
          ic3 <- cluster_sum_ic(cbind(-1/phat^2 * (subset-phat)), id_data)
          al <- align_ic(ics = list(ic1, ic2, ic3),
                         ids = list(uid_data, uid_model, uid_data))
          ic_theta <- (al$ic[[1]] + al$ic[[2]])/phat + rbind(pp)%x%al$ic[[3]]
          pp <- pp/phat
        } else {
          al <- align_ic(list(ic1, ic2), list(uid_data, uid_model))
          ic_theta <- al$ic[[1]] + al$ic[[2]]
        }
        uid <- al$id
        rownames(ic_theta) <- uid
        idstack <- uid
        V <- var_ic(ic_theta)
      }
    }
  }

  if (id_drop) { # id=NULL: remove id (index) and rownames of IF
    idstack <- NULL
    if (!is.null(ic_theta)) rownames(ic_theta) <- NULL
  }

  df_mod <- NULL
  if (inherits(x, "lm") && family(x)$family == "gaussian"
      && !missing(vcov)) {
    # defaults to t-distribution when calculating p-values with model-based SEs
    df_mod <- x$df.residual
  }
  if (is.null(V)) {
    res <- cbind(pp, NA, NA, NA, NA)
  } else {
    if (length(pp)==1)
      res <- rbind(c(pp, diag(V)^0.5))
    else
      res <- cbind(pp, diag(V)^0.5)
  }
  res <- estimate_coefmat(res[, 1], res[, 2], df=df_mod, level=0.95, null=0)
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
                        IC=NULL, print=print, id=idstack, df=df_mod),
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
  res$n <- nrow(data)
  res$ncluster <- if (!is.null(ic_theta)) nrow(ic_theta) else nrow(data)
  res$derivative <- derivative
  res <- structure(res, class="estimate")

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

  if (!is.null(attributes(x)$extra)) {
    cat("\n")
    print(attributes(x)$extra)
  }

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


################################################################################
## Helper functions for estimate.default
################################################################################

## Sum influence function contributions within clusters (first-appearance order
## of 'id') and rescale such that var_ic gives the cluster-robust variance (see
## IC.default). Attributes of 'ic' (e.g. 'bread') are kept, and 'N' (number of
## observations) is set if missing.
cluster_sum_ic <- function(ic, id) {
  if (anyNA(id)) stop("Missing values in 'id'")
  atr <- attributes(ic) # before cbind, which drops attributes
  atr <- atr[setdiff(names(atr), c("dim", "dimnames", "names"))]
  ic <- cbind(ic)
  res <- rowsum(ic, group = id, reorder = FALSE)
  res <- res * NROW(res) / length(id)
  attributes(res)[names(atr)] <- atr
  if (is.null(attr(res, "N"))) attr(res, "N") <- NROW(ic)
  res
}

## Convert an 'id' specification (vector, formula or column name evaluated in
## 'data', or a logical scalar giving 'default') to a vector of length 'n'.
## An 'id' matching the data before removal of missing values (x$na.action) is
## reduced accordingly.
resolve_id <- function(id, data, n, x = NULL, default) {
  if (is.logical(id) && length(id) == 1) return(default)
  if (inherits(id, "formula")) id <- interaction(get_all_vars(id, data))
  if (is.character(id) && length(id) == 1 && n != 1)
    id <- data[, id, drop = TRUE]
  if (length(id) != n) {
    na <- if (is.list(x)) x$na.action
    if (is.null(na) || length(id) != length(na) + n)
      stop("Dimensions of 'id' (", length(id), ") and ",
           "influence function/data (", n, ") does not agree")
    warning("Applying na.action")
    id <- id[-na]
  }
  id
}

## Ids attached to the rows of the influence function 'ic' of 'x': the id of an
## 'estimate' object (keeps the type of the ids), else rownames of 'ic'
ic_ids <- function(x, ic) {
  for (key in list(if (inherits(x, "estimate")) index(x), rownames(ic))) {
    if (length(key) == NROW(ic)) return(key)
  }
  NULL
}

## Align influence functions 'ics' (list) with row ids 'ids' (list of unique
## ids) on the union of ids (first-appearance order). Each IF is zero outside
## its own ids and rescaled by length(union)/length(own ids) (inverse
## probability of observation, see the section "Estimators computed on
## different subsets" in vignette("influencefunction")).
## Returns the list of aligned IFs and the ids of the union.
align_ic <- function(ics, ids) {
  uid <- unique(unlist(ids, use.names = FALSE))
  ics <- Map(function(ic, id) {
    ic <- cbind(ic)
    res <- matrix(0, nrow = length(uid), ncol = ncol(ic),
                  dimnames = list(NULL, colnames(ic)))
    res[match(id, uid), ] <- ic
    res * length(uid) / length(id)
  }, ics, ids)
  list(ic = ics, id = uid)
}

## Ids used when averaging a transformation over 'data' (standardization).
## Returns the ids of the rows of 'data', the model IF aggregated within the
## clusters linked to these ids, and the (unique) model ids.
average_align_ids <- function(x, data, ic, id = NULL) {
  rn <- rownames(data)
  if (is.null(rn)) rn <- as.character(seq_len(NROW(data)))
  id <- average_data_id(x, data, id, rn)
  cl <- average_model_id(x, ic, id, rn)
  list(id_data = id, ic = cluster_sum_ic(ic, cl), id_model = unique(cl))
}

## Ids of the rows of 'data'. Default: index(x) for 'estimate' objects of
## matching length, else rownames (see resolve_id for other specifications).
average_data_id <- function(x, data, id, rn) {
  N <- NROW(data)
  if (is.null(id)) {
    idx <- if (inherits(x, "estimate")) index(x)
    id <- if (length(idx) == N) idx else rn
  }
  id <- resolve_id(id, data, n = N, x = x, default = rn)
  if (is.factor(id)) id <- as.character(id) # avoid integer codes in c(...)
  id
}

## Ids of the rows of the model IF, linked to the ids of 'data' ('id').
## Rules:
##  1. 'estimate' object with ids: used as is (same id space as 'data';
##     non-overlapping ids are treated as independent observations)
##  2. rows of other model objects: linked via the rownames of 'data'
##  3. no ids: positional matching with the rows of 'data'
average_model_id <- function(x, ic, id, rn) {
  key <- ic_ids(x, ic)
  if (inherits(x, "estimate") && !is.null(key)) return(key)
  if (is.null(key)) {
    if (NROW(ic) == length(id)) return(id)
  } else {
    pos <- match(key, rn)
    if (!anyNA(pos)) return(id[pos])
  }
  stop("Unable to link the model influence function to 'data'. ",
       "Supply the ids with 'estimate(x, id=...)'")
}
