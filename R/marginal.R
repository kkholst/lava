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
#' @param f function with the named argument `p` (parameter vector) and
#'   the optional arguments `data` and `object` (the model object). Should
#'   return a vector, matrix or a (named) list of vectors with one value per
#'   row of `data`, or values for some of the rows named by `rownames(data)`
#'   (with `id` given per row of `data`, the average is then conditioned on
#'   these rows, as with `subset`). If
#'   the result has an attribute `"grad"` it is used as the Jacobian of `f`
#'   with respect to `p`.
#' @param data `data.frame` over which `f` is averaged. Defaults to
#'   `model.frame(object)` (which is not available for `estimate` objects).
#' @param id (optional) ids of the values returned by `f` (the g-computation
#'   part). Can be a vector (one value per value of `f`), or refer to the
#'   rows of `data` (a column name, a one-sided formula evaluated in `data`,
#'   or a vector with one value per row of `data`), in which case the values
#'   of `f` are linked to the rows of `data` by their names
#'   (`rownames(data)`). Defaults to the names of the values of `f`, and
#'   otherwise `rownames(data)`. `id = NULL` uses the default ids but removes
#'   the id (index) from the returned object.
#' @param subset (optional) logical vector, expression evaluated in `data`, or
#'   column name. The average is then conditioned on the subpopulation where
#'   `subset` is `TRUE` (conditional marginal estimate). Linked to the values
#'   of `f` via `id`.
#' @param labels (optional) character vector of parameter names
#' @param ... additional arguments passed to `f`
#' @return `estimate` object
#' @details
#'
#' The influence function of the estimate is
#' \deqn{\mathrm{IC}_\Psi(Z; P) = f(X;\theta) - \Psi +
#' [E\nabla_\theta f(bX;\theta)]\,\phi(Z; P)}
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
#' `estimate(object, f, average=TRUE, ...)` is equivalent to
#' `marginal(object, f, ...)`.
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
  if (!is.function(f)) stop("'f' must be a function")
  form <- names(formals(f))
  if (!("p" %in% form)) stop("'f' must have an argument 'p'")
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
    val <- do.call(f, args)
    if (is.list(val)) return(do.call("cbind", val))
    structure(cbind(val, deparse.level = 0), grad = attr(val, "grad"))
  }
  val <- fval(pp)
  D <- attr(val, "grad")
  if (is.null(D)) D <- numDeriv::jacobian(fval, pp)
  N <- NROW(val)
  if (N <= 1L)
    stop("'f' must return one value per observation (row of 'data')")
  nn <- colnames(val)
  if (is.null(nn)) nn <- paste0("p", seq_len(NCOL(val)))

  ## Ids of the values of 'f': 'id', else the names of the values, else
  ## rownames(data). 'id' (and 'subset') may also refer to the rows of 'data'
  ## (column name, or one value per row of 'data'), linked to the values by
  ## their names (rownames(data)). Rows without a value of 'f' are then
  ## excluded from the average (conditional average).
  nd <- NROW(data)
  pos <- if (!is.null(rownames(val)) && !is.null(rownames(data)))
           match(rownames(val), rownames(data))
  if (anyNA(pos)) pos <- NULL
  id_drop <- !missing(id) && is.null(id)
  expanded <- FALSE
  if (missing(id) || is.null(id)) {
    id <- rownames(val)
    if (is.null(id) && !is.null(data)) {
      if (nd != N)
        stop("Unable to link the values of 'f' (", N, ") to the rows of ",
             "'data' (", nd, "). Supply 'id' or name the values of 'f'")
      id <- rownames(data)
    }
    if (is.null(id)) id <- seq_len(N)
  } else {
    per_row <- inherits(id, "formula") || (is.character(id) && length(id) == 1)
    if (per_row) id <- resolve_id(id, data, n = nd, x = object)
    per_row <- per_row || (length(id) == nd && length(id) != N)
    if (per_row && !is.null(pos)) { # expand values to the rows of 'data'
      val <- expand_rows(val, pos, nd)
      D <- expand_rows(D, as.vector(outer(pos, (seq_len(ncol(val)) - 1) * nd,
                                          "+")), nd * ncol(val))
      N <- nd
      expanded <- TRUE
    } else if (length(id) != N) {
      id <- resolve_id(id, data, n = N, x = object) # na.action / error
      pos <- NULL
    }
  }
  if (is.factor(id)) id <- as.character(id) # avoid integer codes in c(...)
  obs <- if (expanded) seq_len(N) %in% pos else rep(TRUE, N)

  if (!missing(subset)) {
    s <- resolve_subset(substitute(subset), data, parent.frame())
    if (!expanded && !is.null(pos) && length(s) == nd) s <- s[pos]
    if (length(s) != N)
      stop("Length of 'subset' (", length(s), ") does not agree with the ",
           "number of values of 'f' (", N, ")")
    obs <- obs & s
  }

  res <- average_ic(val, D, ic, id_data = id,
                    id_model = average_model_id(e, ic),
                    subset = if (!all(obs)) obs)
  names(res$coef) <- nn
  colnames(res$ic) <- nn
  out <- estimate(coef = res$coef, IC = res$ic, id = res$id)
  if (!is.null(labels)) out <- labels(out, labels)
  if (id_drop) index(out) <- NULL
  out$n <- sum(obs)
  out$derivative <- attr(val, "grad")
  out
}

## Expand the rows of matrix 'x' to 'n' rows (rows 'pos'; zero elsewhere)
expand_rows <- function(x, pos, n) {
  res <- matrix(0, n, NCOL(x))
  res[pos, ] <- x
  res
}

## Evaluate 'subset' (an expression, evaluated in 'env' and otherwise in
## 'data', or a column name of 'data') to a logical vector
resolve_subset <- function(expr, data, env) {
  s <- tryCatch(eval(expr, envir = env), error = function(err) {
    if (is.null(data)) stop(err)
    eval(expr, envir = data, enclos = env)
  })
  if (is.character(s) && length(s) == 1) s <- data[, s, drop = TRUE]
  if (is.numeric(s)) s <- s > 0
  s
}


## Ids attached to the rows of the influence function 'ic' of 'x': the id of an
## 'estimate' object (keeps the type of the ids), else rownames of 'ic'
ic_ids <- function(x, ic) {
  for (key in list(if (inherits(x, "estimate")) index(x), rownames(ic))) {
    if (length(key) == NROW(ic)) return(key)
  }
  NULL
}

## Ids of the rows of the model IF: index(x) for 'estimate' objects, else the
## rownames of the IF (for model objects: the rownames of the model frame).
## The ids are in the same id space as the ids of 'data'; non-overlapping ids
## are treated as independent observations.
average_model_id <- function(x, ic) {
  key <- ic_ids(x, ic)
  if (!is.null(key)) return(key)
  stop("Unable to link the model influence function to 'data'. ",
       "Supply the ids with 'estimate(x, id=...)'")
}

## Influence function of the empirical average of f(X; theta).
## val: N x k matrix of f(X_i; theta); D: Jacobian (N*k x p, column-major in
## val) of val wrt theta; ic: row-level model IF; id_data, id_model: ids of
## the rows of 'val' and 'ic'; subset: logical vector (conditional average).
## The empirical term (on the ids of the data) and the model term (on the ids
## of the model) are aligned on the union of the ids (see 'align_ic').
## Returns the estimates, the IF and the ids of the rows of the IF.
average_ic <- function(val, D, ic, id_data, id_model, subset = NULL) {
  val <- cbind(val)
  N <- nrow(val)
  k <- ncol(val)
  w <- if (is.null(subset)) rep(1, N) else as.numeric(subset)
  val <- val * w
  ## Average derivative wrt theta (k x p)
  D <- rbind(D)
  D0 <- matrix(nrow = k, ncol = ncol(D))
  for (i in seq_len(k)) {
    D0[i, ] <- colMeans(D[seq_len(N) + (i - 1) * N, , drop = FALSE] * w)
  }
  pp <- colMeans(val)
  ## Model term, D phi, on the ids of the model
  ic2 <- cluster_sum_ic(ic, id_model) %*% t(D0)
  uid_model <- unique(id_model)
  ## Empirical averaging term, f(X; theta) - Psi, on the ids of 'data'
  ic1 <- cluster_sum_ic(val - matrix(pp, N, k, byrow = TRUE), id_data)
  uid_data <- unique(id_data)
  if (!any(uid_model %in% uid_data)) {
    message("Assuming independence between model iid decomposition and new data frame") #nolint
  }
  if (!is.null(subset)) { ## Conditional estimate
    phat <- mean(w)
    ic3 <- cluster_sum_ic(cbind(-1 / phat^2 * (w - phat)), id_data)
    al <- align_ic(list(ic1, ic2, ic3), list(uid_data, uid_model, uid_data))
    ic <- (al$ic[[1]] + al$ic[[2]]) / phat + al$ic[[3]] %*% rbind(pp)
    pp <- pp / phat
  } else {
    al <- align_ic(list(ic1, ic2), list(uid_data, uid_model))
    ic <- al$ic[[1]] + al$ic[[2]]
  }
  rownames(ic) <- NULL
  list(coef = pp, ic = ic, id = al$id)
}



## Align influence functions 'ics' (list) with row ids 'ids' (list of unique
## ids) on the union of ids (first-appearance order). Each IF is zero outside
## its own ids and rescaled by length(union)/length(own ids) (inverse
## probability of observation, see the section "Estimators computed on
## different subsets" in vignette("influencefunction")).
## Returns the list of aligned IFs and the ids of the union.
align_ic <- function(ics, ids) {
  ## union of ids. Starting from the first set keeps its type (e.g. numeric
  ## ids) when the other sets are contained in it
  uid <- unique(ids[[1]])
  for (i in ids[-1]) {
    i <- unique(i)
    new <- is.na(match(i, uid))
    if (any(new)) uid <- c(uid, i[new])
  }
  ics <- Map(function(ic, id) {
    ic <- cbind(ic)
    res <- matrix(0, nrow = length(uid), ncol = ncol(ic),
                  dimnames = list(NULL, colnames(ic)))
    res[match(id, uid), ] <- ic
    res * length(uid) / length(id)
  }, ics, ids)
  list(ic = ics, id = uid)
}
