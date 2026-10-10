#' Marginal (standardized) estimates via g-computation
#'
#' Estimates the marginal (standardized) parameter
#' \eqn{\Psi = E\{f(X; \theta)\}} by the empirical average of `f(p, data)`
#' over the rows of `data` (g-computation), where the parameter \eqn{\theta}
#' is estimated by `object`. The influence function of the estimate
#' accounts for both the empirical averaging and the uncertainty of
#' \eqn{\widehat\theta}.
#'
#' @param object model object (`glm`, `lvmfit`, ...) or `estimate` object
#'   providing the parameter estimates and their influence function.
#' @param f function where the first argument is the parameter vector and
#'   the optional named arguments `data` and `object` (the model object, e.g.
#'   an `estimate` object. Should return a vector, matrix or a (named) list of
#'   vectors with one value per observation (row of `data`). If the result has
#'   an attribute `"grad"` it is used as the Jacobian of `f` with respect to `p`
#' @param data `data.frame` over which `f` is averaged. Defaults to
#'   `model.frame(object)` (which is not available for `estimate` objects).
#' @param id (optional) ids of the values returned by `f` (the g-computation
#'   part): a vector with one value per value of `f`, or a column name or
#'   one-sided formula evaluated in `data` (requiring one value of `f` per
#'   row of `data`). Defaults to the names of the values of `f`, and
#'   otherwise `rownames(data)`. `id = NULL` uses the default ids but removes
#'   the id (index) from the returned object.
#' @param subset (optional) logical vector (one value per value of `f`),
#'   expression evaluated in `data` (columns of `data` take precedence over
#'   variables in the calling environment), or column name. The average is
#'   then conditioned on the subpopulation where `subset` is `TRUE`
#'   (conditional marginal estimate).
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
#' g <- glm(y ~ z + x, data = d, family = binomial())
#'
#' ## Standardization (g-computation)
#' f <- function(p, data)
#'   list(p0 = expit(p["(Intercept)"] + p["z"] * data[, "z"]),
#'        p1 = expit(p["(Intercept)"] + p["x"] + p["z"] * data[, "z"]))
#' e <- marginal(g, f)
#' e
#' estimate(e, diff)
#'
#' ## Model estimated on a subset, standardized over the full data
#' d$id <- paste0("i", seq_len(nrow(d)))
#' d$w <- rbinom(nrow(d), 1, 0.5)
#' d1 <- subset(d, w == 1)
#' g1 <- glm(y ~ x + z, data = d1, family = binomial)
#' e1 <- estimate(g1, id = d1$id)
#' marginal(e1, f, data = d, id = "id")
#'
#' ## Conditional marginal effects (subset) and clusters
#' d$cl <- rep(seq(nrow(d) / 4), each = 4)
#' marginal(estimate(g, id = d$cl),
#'          function(p, data) expit(p[1] + p["z"] * data[, "z"]),
#'          data = d, subset = z > 0, id = "cl")
marginal <- function(object, f, data, id, subset, labels = NULL, ...) {
  form <- checkarg(object, f, "p")
  ## Model: parameter estimates and (row-level) influence function
  e <- if (inherits(object, "estimate")) object else estimate(object)
  ic <- IC(e)
  if (is.null(ic)) stop("Influence function of 'object' needed")
  pp <- coef(e)
  if (missing(data))
    data <- tryCatch(model.frame(object), error = function(...) NULL)

  ## Values of 'f' (N x k matrix) and Jacobian (N*k x p)
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
  D <- grad <- attr(val, "grad")
  if (is.null(D)) D <- numDeriv::jacobian(fval, pp,
                                          method=lava.options()$Dmethod)
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
    if (is.null(id)) id <- seq_len(N)
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
    subset <- resolve_subset(substitute(subset), data, parent.frame())
    if (length(subset) != N)
      stop("Length of 'subset' (", length(subset), ") does not agree with ",
           "the number of values of 'f' (", N, ")")
  } else {
    subset <- NULL
  }

  colnames(val) <- nn
  out <- average_estimate(val, D, ic, id_data = id, id_model = id_model,
                          subset = subset)
  out <- labels(out, if (is.null(labels)) nn else labels)
  if (id_drop) index(out) <- NULL
  out$n <- if (is.null(subset)) N else sum(subset)
  out$derivative <- grad
  out
}

# Evaluate 'subset' (an expression evaluated in 'data' (enclosed by 'env'),
# or a column name of 'data') to a logical vector. Columns of 'data' take
# precedence over variables in 'env' (as in 'subset', 'with')
resolve_subset <- function(expr, data, env) {
  s <- if (!is.list(data)) eval(expr, envir = env) else
         eval(expr, envir = data, enclos = env)
  if (is.character(s) && length(s) == 1) s <- data[, s, drop = TRUE]
  if (is.numeric(s)) s <- s > 0
  s
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
average_estimate <- function(val, D, ic, id_data, id_model, subset = NULL) {
  val <- cbind(val)
  N <- nrow(val)
  k <- ncol(val)
  w <- if (is.null(subset)) rep(1, N) else as.numeric(subset)
  psi <- colSums(val * w) / sum(w)
  G <- rowsum(rbind(D) * w, rep(seq_len(k), each = N), reorder = FALSE) /
    sum(w)
  if (!any(unique(id_model) %in% unique(id_data))) {
    message("Assuming independence between model iid decomposition and new data frame") #nolint
  }
  e_data <- estimate(coef = psi, id = id_data,
                     IC = w * (val - matrix(psi, N, k, byrow = TRUE)) / mean(w))
  e_model <- estimate(coef = 0 * psi, IC = ic %*% t(G), id = id_model)
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
