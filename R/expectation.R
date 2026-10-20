#' Mean of (estimated) function via g-computation
#'
#' Estimates the mean parameter
#' \eqn{\Psi = E{f(X; \theta)}} by the empirical average of f(p, data)
#' over the rows of data, where the parameter \eqn{\theta}
#' is estimated by object. The influence function of the estimate
#' accounts for both the empirical averaging and the uncertainty of
#' \eqn{\widehat\theta}. Standardization and g-computation are special
#' applications of this function.
#'
#' @param object model object (`glm`, `lvmfit`, ...) or `estimate` object
#'   providing the parameter estimates and their influence function.
#'   Alternatively, an evaluated expression of a function with values in `data`
#'   can be provided to calculate the simple mean (see examples).
#' @param f function with named argument `p` providing the parameter vector and
#'   the optional named arguments `data` and `object` (the model object, e.g. an
#'   `estimate` object). Should return a vector, matrix or a (named) list of
#'   vectors with one value per observation (row of `data` unless `data` is
#'   missing). If the result has an attribute `"grad"` it is used as the
#'   Jacobian of `f` with respect to `p`
#' @param data `data.frame` over which `f` is averaged. Defaults to
#'   `model.frame(object)` (which is not available for `estimate` objects).
#' @param id (optional) ids of the values returned by `f` (the g-computation
#'   part): a vector with one value per value of `f`, or a column name or
#'   one-sided formula evaluated in `data` (requiring one value of `f` per row
#'   of `data`). Defaults to the names of the values of `f`, and otherwise
#'   `rownames(data)`. `id = NULL` uses the default ids but removes the id
#'   (index) from the returned object.
#' @param subset (optional) subpopulation on which the average is conditioned
#'   (conditional marginal mean estimate). Either a logical vector (or numeric
#'   0/1 vector) with one value per value of `f`, a column name of `data`, or a
#'   one-sided formula evaluated in `data` (e.g., `~ a == 1 & z > 0`). Note
#'   that `subset` uses standard evaluation, i.e., `subset = a == 1` is
#'   evaluated in the calling environment and not in `data`. As for model
#'   formulas (see [stats::model.frame]), the variables of a formula `subset`
#'   are looked up in the columns of `data` and then in the environment where
#'   the formula was created (not the environment from which `expectation` is
#'   called).
#' @param labels (optional) character vector of parameter names
#' @param ... additional arguments passed to `f`
#' @return `estimate` object
#' @details
#'
#' The influence function of the estimate is
#' \deqn{\mathrm{IC}_\Psi(Z; P) = f(X;\theta) - \Psi +
#' [E\nabla_\theta f(X;\theta)]\,\phi(Z; P)}
#' where \eqn{\phi} is the influence function of \eqn{\widehat\theta}.
#'
#' The two terms are identified by different ids: the first term by `id` (the
#' rows of `data`), and the second term by the ids of the model, i.e.,
#' `index(object)` for `estimate` objects, or the row names of the influence
#' function (for model objects such as `glm` the row names of the model
#' frame). The model may therefore be estimated on another (e.g., smaller or
#' partly overlapping) dataset than `data`. The two terms are aligned on the
#' union of ids: each term is zero for ids outside its own support and
#' rescaled by the inverse proportion of observed ids (as in
#' [merge.estimate]). If there are no common ids, the model estimate and
#' `data` are treated as independent. Use `estimate(object, id=...)` to
#' assign ids (or clusters) to the model.
#'
#' @seealso [estimate.default], `vignette("influencefunction", package =
#'   "lava")`
#' @export
#' @examples
#' m <- lvm(y ~ x + z)
#' distribution(m, y ~ x) <- dist_bernoulli("logit")
#' d <- sim(m, 1000, seed = 1)
#' fit <- glm(y ~ z + x, data = d, family = binomial())
#'
#' ## Standardization (g-computation)
#' f <- function(p, data)
#'   list(p0 = expit(p["(Intercept)"] + p["z"] * data[, "z"]),
#'        p1 = expit(p["(Intercept)"] + p["x"] + p["z"] * data[, "z"]))
#' e <- expectation(f, fit)
#' e
#' estimate(e, diff)
#'
#' ## Model estimated on a subset, standardized over the full data
#' d$id <- paste0("i", seq_len(nrow(d)))
#' d$w <- rbinom(nrow(d), 1, 0.5)
#' d1 <- subset(d, w == 1)
#' g1 <- glm(y ~ x + z, data = d1, family = binomial)
#' e1 <- estimate(g1, id = d1$id)
#' expectation(f, e1, data = d, id = "id")
#'
#' ## Conditional marginal mean effects (subset) and clusters
#' d$cl <- rep(seq(nrow(d) / 4), each = 4)
#' expectation(function(p, data) expit(p[1] + p["z"] * data[, "z"]),
#'             estimate(fit, id = d$cl),
#'             data = d, subset = ~ z > 0, id = "cl")
#'
#' ## Simple mean of known function
#' expectation(~ 1 / (1+z^2), data=d, subset = ~x, labels="p")
expectation <- function(f, object, data, id, subset, labels = NULL, ...) {

  if (!missing(subset)) {
    subset <- tryCatch(subset, error = function(e) {
      stop("Unable to evaluate 'subset': ", conditionMessage(e),
           ". To evaluate an expression in 'data' use a one-sided ",
           "formula, e.g., `subset = ~ a == 1`", call. = FALSE)
    })
  }

  if (missing(object)) {
    if (missing(data)) data <- NULL
    out <- expect(
      expr = f,
      data = data,
      envir = parent.frame(),
      id = id,
      quoted = TRUE,
      subset = subset,
      ...
    )
    if (!is.null(labels))
      out <- labels(out, labels)
    return(out)
  }

  # Otherwise we calculate E[f(Z; theta)], with theta estimated by object
  form <- checkarg(object, f, "p")
  ## Model: parameter estimates and (row-level) influence function
  e <- if (inherits(object, "estimate")) object else estimate(object)
  ic <- IC(e)
  if (is.null(ic)) stop("Influence function of 'object' needed")
  pp <- coef(e)
  if (missing(data))
    data <- tryCatch(model.frame(object), error = function(...) NULL)

  ## Values of 'f' (N x k matrix)
  dots <- list(...)
  fval <- function(p) {
    args <- c(list(p = p, data = data, object = object), dots)
    if (!("..." %in% form)) args <- args[names(args) %in% form]
    ## args <- c(list(p), args)
    val <- do.call(f, args)
    if (is.list(val)) return(do.call("cbind", val))
    structure(cbind(val, deparse.level = 0), grad = attr(val, "grad"))
  }
  val <- fval(pp)

  N <- NROW(val)
  if (N <= 1L)
    stop("'f' must return one value per observation (row of 'data')")
  nn <- colnames(val)
  if (is.null(nn)) nn <- paste0("p", seq_len(NCOL(val)))

  ## Ids of the values of 'f': 'id', else the names of the values, else
  ## rownames(data)
  id_drop <- !missing(id) && is.null(id)
  if (missing(id) || is.null(id)) {
    id <- rownames(val)
    if (is.null(id) && !is.null(data)) {
      if (NROW(data) != N)
        stop("Unable to link the values of 'f' (", N, ") to the rows of ",
             "'data' (", NROW(data), "). Supply 'id' or name the values of 'f'")
      id <- rownames(data)
    }
    ## if (is.null(id)) id <- seq_len(N)
  } else {
    id <- resolve_id(id, data, n = N, x = object)
  }
  if (is.factor(id)) id <- as.character(id) # avoid integer codes in c(...)

  ## Ids of the rows of the model IF: index(x) for 'estimate' objects, else
  ## the rownames of the IF (for model objects: rownames of the model frame)
  id_model <- ic_ids(e, ic)
  if (is.null(id_model))
    stop("Unable to link the model influence function to 'data'. ",
         "Supply the ids with 'estimate(x, id=...)'")

  if (!missing(subset)) {
    subset <- resolve_subset(subset, data, N)
  } else {
    subset <- NULL
  }

  w <- if (is.null(subset)) rep(1, N) else as.numeric(subset)
  D <- attr(val, "grad") # optional analytic Jacobian (N*k x p)
  val <- cbind(val)
  k <- NCOL(val)
  colnames(val) <- nn
  ## Derivative of the (weighted) average, k x p (dim(f) x dim(theta))
  G <- if (!is.null(D)) {
    rowsum(rbind(D) * w, rep(seq_len(k), each = N), reorder = FALSE) / sum(w)
       } else {
         meanf <- function(theta) colSums(fval(theta) * w) / sum(w)
         rbind(numDeriv::jacobian(meanf,
                             pp, method = lava.options()$Dmethod))
  }
  out <- average_estimate(val, G, ic,
                          id_data = id, id_model = id_model,
                          w = w)
  out <- labels(out, if (is.null(labels)) nn else labels)
  if (id_drop) index(out) <- NULL
  out$n <- if (is.null(subset)) N else sum(subset)
  out
}

# Calculate E[f(Z)]
expect <- function(expr, data = NULL, id, subset, ...) {
  expr <- tryCatch(expr, error = function(e) {
      stop("Unable to evaluate 'expr': ", conditionMessage(e),
           ". To evaluate an expression in 'data' use a one-sided ",
           "formula, e.g., `expr = ~ y`", call. = FALSE)
  })
  if (inherits(expr, "formula")) {
    if (length(expr) != 2L)
      stop("'expr' formula must be one-sided, e.g., `~ y`")
    expr <- eval(expr[[2L]], envir = data, enclos = environment(expr))
  } else if (is.character(expr) && length(expr) == 1L) {
    if (is.null(data) || !(expr %in% names(data)))
      stop("'expr' column '", expr, "' not found in 'data'")
    expr <- data[[expr]]
  }
  z <- cbind(expr, deparse.level = 0)
  N <- NROW(z)
  if (!missing(subset)) {
    subset <- resolve_subset(subset, data, N)
  } else {
    subset <- rep(1, N)
  }
  id_drop <- !missing(id) && is.null(id)
  if (missing(id) || is.null(id)) {
    id <- rownames(z)
    if (is.null(id) && !is.null(data)) {
      if (NROW(data) != N)
        stop("Unable to link the values of 'f' (", N, ") to the rows of ",
             "'data' (", NROW(data), "). Supply 'id' or name the values of 'f'")
      id <- rownames(data)
    }
  } else {
    id <- resolve_id(id, data, n = N, x = NULL)
  }
  if (is.factor(id)) id <- as.character(id)
  w <- subset / mean(subset)
  est <- colMeans(w * z)
  ic <- matrix(ncol=ncol(z), nrow=nrow(z))
  for (i in seq_len(ncol(z))) {
    ic[, i] <- w * (z[, i] - est[i])
  }
  res <- estimate(coef=est, IC=ic, id = id, ...)
  if (id_drop) index(res) <- NULL
  res
}

# Resolve the 'subset' argument to a logical vector of length 'n':
#  - one-sided formula: evaluated in 'data' enclosed by environment(formula)
#    (as in model.frame), i.e., columns of 'data' mask variables of the same
#    name in the environment where the formula was created
#  - character string: column of 'data'
#  - logical vector, or numeric vector with values in {0, 1}
resolve_subset <- function(subset, data, n) {
  if (is.null(subset)) return(NULL)
  if (!is.list(data)) data <- NULL
  if (inherits(subset, "formula")) {
    if (length(subset) != 2L)
      stop("'subset' formula must be one-sided, e.g., `~ a == 1`")
    subset <- eval(subset[[2L]], envir = data, enclos = environment(subset))
  } else if (is.character(subset) && length(subset) == 1L) {
    if (is.null(data) || !(subset %in% names(data)))
      stop("'subset' column '", subset, "' not found in 'data'")
    subset <- data[[subset]]
  }
  if (is.numeric(subset) && !anyNA(subset)) {
    if (!all(subset %in% c(0, 1)))
      stop("Numeric 'subset' must only contain the values 0 and 1")
    subset <- subset == 1
  }
  if (!is.logical(subset))
    stop("'subset' must be a logical vector, a column name of 'data', or a ",
         "one-sided formula evaluated in 'data' (e.g., `~ a == 1`)")
  if (anyNA(subset)) stop("Missing values in 'subset'")
  if (length(subset) != n)
    stop("Length of 'subset' (", length(subset), ") does not agree with ",
         "the number of values of 'f' (", n, ")")
  subset
}

# Ids attached to the rows of the influence function 'ic' of 'x': the id of an
# 'estimate' object (keeps the type of the ids), else rownames of 'ic'
ic_ids <- function(x, ic) {
  for (key in list(if (inherits(x, "estimate")) index(x), rownames(ic))) {
    if (length(key) == NROW(ic)) return(key)
  }
  NULL
}

# Estimate of the average E[f(X; theta)] (or conditional_E[f(X; theta)|W=1])
average_estimate <- function(val, expectdf, ic,
                             id_data, id_model, w = NULL) {
  val <- cbind(val)
  N <- nrow(val)
  k <- ncol(val)
  if (is.null(w)) w <- rep(1, N)
  psi <- colSums(val * w) / sum(w)
  if (!any(unique(id_model) %in% unique(id_data))) {
    message("Assuming independence between model iid decomposition and new data frame") #nolint
  }
  e_data <- estimate(coef = psi, id = id_data,
                     IC = w * (val - matrix(psi, N, k, byrow = TRUE)) / mean(w))
  e_model <- estimate(coef = 0 * psi, IC = ic %*% t(expectdf), id = id_model)
  e_data + e_model # IC automatically aligned via merge.estimate
}

checkarg <- function(object, f, arg=NULL) {
  if (!is.function(f)) stop("'f' must be a function")
  generic <- utils::isS3stdGeneric(f)
  if (isTRUE(generic)) {
    method <- utils::getS3method(names(generic), class(object)[1L])
    form <- names(formals(method))
  } else {
    form <- names(formals(f))
  }
  if (!is.null(arg) && !(arg %in% form)) {
    stop("'f' must have an argument '", arg, "'")
  }
  return(form)
}
