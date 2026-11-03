#' Expand long-format data to all subject and visit combinations
#'
#' Returns a long-format data.frame with exactly one row per subject and visit
#' (`levels`), sorted by subject (in order of first appearance) and visit.
#' Rows that are not present in `data` are added by copying all columns from
#' the previous existing row of the subject (last observation carried
#' forward), or from the next existing row if there is no previous row. The
#' response of added rows is set to `NA`.
#' @title Expand long-format data
#' @param data (data.frame) data in long format
#' @param id (character) name of the subject variable
#' @param time (character) name of the visit variable. In the returned data
#'   this variable is a factor with levels `levels`.
#' @param response (character) optional name of the response variable which is
#'   set to `NA` in added rows
#' @param levels (vector) visit levels. Defaults to the levels of the visit
#'   variable if it is a factor, and otherwise its sorted unique values.
#' @return data.frame with `length(levels)` rows per subject
#' @export
#' @examples
#' d <- data.frame(id = c(1, 1, 2), visit = c(1, 2, 2), y = 1:3, x = 4:6)
#' complete_visits(d, "id", "visit", "y")
complete_visits <- function(data, id, time, response = NULL, levels = NULL) {
  if (is.null(levels)) {
    levels <- if (is.factor(data[[time]])) {
      levels(data[[time]])
    } else {
      sort(unique(data[[time]]))
    }
  }
  ids <- unique(as.character(data[[id]]))
  i <- match(as.character(data[[id]]), ids)
  j <- match(as.character(data[[time]]), as.character(levels))
  if (anyNA(j)) stop("'data' contains visits not in 'levels'.")
  if (anyDuplicated(cbind(i, j))) stop("Duplicated subject/visit rows.")
  idx <- matrix(NA_integer_, length(ids), length(levels))
  idx[cbind(i, j)] <- seq_len(nrow(data))
  src <- idx
  for (k in seq_along(levels)[-1]) { # LOCF
    src[, k] <- ifelse(is.na(src[, k]), src[, k - 1], src[, k])
  }
  for (k in rev(seq_along(levels))[-1]) { # first visit(s) missing: next row
    src[, k] <- ifelse(is.na(src[, k]), src[, k + 1], src[, k])
  }
  res <- data[as.vector(t(src)), , drop = FALSE]
  res[[time]] <- factor(rep(levels, length(ids)), levels = levels)
  if (!is.null(response)) res[[response]][is.na(as.vector(t(idx)))] <- NA
  rownames(res) <- NULL
  return(res)
}

#' Predictions at a single visit from an mmrm model
#'
#' For each subject in `newdata`, returns the outcome at visit `time` if it is
#' observed, and otherwise the conditional mean given the observed outcomes of
#' the subject. With \eqn{O} the observed visits of subject \eqn{i},
#' \deqn{E(Y_{it} \mid Y_{iO}) = x_{it}^\top\beta +
#' \Sigma_{tO}\Sigma_{OO}^{-1}(Y_{iO} - X_{iO}\beta),}
#' where \eqn{\Sigma} is the (group-specific) covariance matrix of all visits
#' (the marginal mean \eqn{x_{it}^\top\beta} if no outcomes are observed).
#' This corresponds to `predict(fit, newdata, conditional = TRUE)` of the mmrm
#' package.
#'
#' `newdata` is first expanded to all subject and visit combinations with
#' [complete_visits()].
#' @title Predictions at a single visit from mmrm models
#' @param object (estimate.mmrm) estimate of an mmrm model obtained with
#'   `estimate(fit, sigma = TRUE)`.
#' @param newdata (data.frame) long-format data with the subject and visit
#'   variables, covariates and (optionally) the response. Defaults to the data
#'   used to fit the model.
#' @param time (integer) index of the visit to predict (default last visit).
#' @param p (numeric) parameter vector on the scale of `coef(object)`: mean
#'   parameters followed by the upper-triangular elements (column-major,
#'   including the diagonal) of the covariance matrix of each group.
#' @param ... additional arguments (not used)
#' @return Named numeric vector (subject IDs) with attributes `observed`
#'   (logical) and `time`.
#' @export
#' @examples
#' if (requireNamespace("mmrm", quietly = TRUE)) {
#'   data(fev_data, package = "mmrm")
#'   fit <- mmrm::mmrm(
#'     FEV1 ~ ARMCD * AVISIT + us(AVISIT | USUBJID),
#'     data = fev_data, reml = FALSE,
#'     optimizer = "nlminb"
#'   )
#'   e <- estimate(fit, sigma = TRUE)
#'   predict(e, newdata = subset(fev_data, USUBJID %in% c("PT1", "PT2")))
#' }
predict.estimate.mmrm <- function(object, newdata = NULL, time = NULL,
                                  p = coef(object), ...) {
  if (!isTRUE(object$sigma)) {
    stop("'object' must be obtained with estimate(fit, sigma = TRUE).")
  }
  stopifnot(length(p) == length(coef(object)))
  fit <- object$fit
  fp <- fit$formula_parts
  object$coef <- p
  cc <- coef(object, list = TRUE) # mean (beta) and covariance (sigma) params
  visits <- levels(fit$tmb_data$full_frame[[fp$visit_var]])
  groups <- "sigma"
  if (fit$tmb_data$n_groups > 1L) groups <- levels(fit$tmb_data$subject_groups)
  Sig <- vec2sigma(cc$sigma, groups = groups, visits = visits)
  K <- length(visits)
  if (is.null(time)) time <- K
  stopifnot(length(time) == 1L, time %in% seq_len(K))
  if (is.null(newdata)) newdata <- fit$tmb_data$data
  resp <- all.vars(fp$model_formula[[2]])
  if (!all(resp %in% names(newdata))) newdata[resp] <- NA_real_
  d <- complete_visits(newdata, fp$subject_var, fp$visit_var, resp[1], visits)

  mu <- matrix(.mmrm_design(fit, d) %*% cc$beta, ncol = K, byrow = TRUE)
  y <- matrix(eval(fp$model_formula[[2]], d), ncol = K, byrow = TRUE)
  r <- y - mu # NA if outcome or covariates are missing
  first <- seq(1, nrow(d), by = K)
  g <- rep(1L, length(first)) # covariance group of each subject
  if (length(Sig) > 1L) g <- as.character(d[[fp$group_var]][first])

  pred <- ifelse(is.na(r[, time]), NA_real_, y[, time])
  cache <- new.env(hash = TRUE, parent = emptyenv())
  observed <- !is.na(pred)
  for (i in which(!observed)) {
    O <- which(!is.na(r[i, ]))
    pred[i] <- mu[i, time]
    if (length(O) == 0L) next
    # hash/key: variance group, missing pattern
    key <- paste0(g[i], ":", paste(O, collapse = ","))
    w <- get0(key, envir = cache, inherits = FALSE) # NULL if key not found
    if (is.null(w)) { # calculate and cache
      S <- Sig[[g[i]]]
      w <- solve(S[O, O, drop = FALSE], S[O, time])
      assign(key, w, envir = cache)
    }
    pred[i] <- pred[i] + sum(w * r[i, O])
  }
  structure(pred, names = as.character(d[[fp$subject_var]][first]),
            observed = observed, time = time)
}

#' @export
coef.estimate.mmrm <- function(object,
                               list = FALSE, ...) {
  p <- object$coef
  if (!list) return(p)
  np <- length(object$fit$beta_est)
  base::list(beta=p[seq_len(np)],
             sigma=p[-seq_len(np)])
}

## Fixed-effects design matrix of an mmrm fit evaluated in 'data' (rows with
## missing predictors are kept as NA rows).
.mmrm_design <- function(fit, data) {
  tt <- stats::delete.response(stats::terms(fit$formula_parts$model_formula))
  mf <- stats::model.frame(tt, data,
                           xlev = mmrm::component(fit, "xlev"),
                           na.action = stats::na.pass)
  X <- stats::model.matrix(tt, mf,
                           contrasts.arg = mmrm::component(fit, "contrasts"))
  return(X[, colnames(mmrm::component(fit, "x_matrix")), drop = FALSE])
}
