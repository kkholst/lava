#' @export
score.mmrm <- function(x,
                       p = NULL,
                       indiv = TRUE,
                       which = c("beta", "theta"),
                       ...) {
  nbeta <- length(coef(x))
  if (is.null(p)) {
    p <- pars(x)
    beta <- p[seq_len(nbeta)]
    theta <- p[-seq_len(nbeta)]
  } else {
    if ("beta" %in% which) {
      beta <- p[seq_len(nbeta)]
      p <- p[-seq_along(beta)]
    }
    if ("theta" %in% which) {
      theta <- p
    }
  }
  U1 <- U2 <- NULL
  if ("beta" %in% which)  U1 <- .mmrm_score_beta(x, beta=beta)
  if ("theta" %in% which) U2 <- .mmrm_score_theta(x, theta=theta)
  res <- cbind(U1, U2)
  ## Attach subject IDs as row names
  rownames(res) <- vapply(.mmrm_subjects(x),
                          `[[`, character(1), "id")
 if (!indiv) return(colSums(res))
  return(res)
}

#' @export
pars.mmrm <- function(x, which = c("beta", "theta"),
                      list = FALSE, p = NULL, ...) {
  beta <- x$beta_est
  theta <- x$theta_est
  if (is.null(names(theta))) {
    names(theta) <- paste0("theta", seq_along(theta))
  }
  if (!is.null(p)) {
    pos <- 0
    if ("beta" %in% which) {
      beta <- structure(p[seq_along(x$beta_est)], names=names(beta))
      pos <- length(beta)
    }
    if ("theta" %in% which) {
      theta <- structure(p[seq_along(x$theta_est) + pos], names=names(theta))
    }
  }
  res <- c()
  if ("beta" %in% which) {
    res <- beta
    if (list) {
      res <- base::list(beta)
      names(res) <- "beta"
    }
  }
  if ("theta" %in% which) {
    if (list) {
      item <- base::list(theta)
      names(item) <- "theta"
      res <- c(res, item)
    } else {
      res <- c(res, theta)
    }
  }
  return(res)
}

#' @export
IC.mmrm <- function(x, ..., numeric=FALSE) {
  pp <- pars(x, ...)
  if (numeric) {
    I <- -numDeriv::jacobian(function(p) {
      score(x, p = p, indiv = FALSE)#, ...)
    }, pp, method = lava::lava.options()$Dmethod)
  } else {
  ## The hessian is block-diagonal and can be extracted
  ## directly from the mmrm / TMB object
    I1 <- solve(x$beta_vcov)
    I2 <- solve(x$theta_vcov) # should equal x$tmb_object$he()
    nn <- paste0("theta", seq_len(nrow(I2)))
    dimnames(I2) <- list(nn, nn)
    I <- lava::blockdiag(I1, I2)
  }
  U <- score(x, indiv=TRUE, ...)
  bread <- Inverse(I) * NROW(U)
  res <- U %*% bread
  colnames(res) <- names(pp)
  return(res)
}

#' @export
estimate.mmrm <- function(x,
                          which = c("beta", "theta"),
                          sigma = FALSE,
                          id = NULL,
                          ...) {
  ic <- IC(x, which = which)
  if (is.null(id)) id <- rownames(ic)
  res <- lava::estimate(coef = pars(x, which = which),
                        IC = ic, id = id)
  if (sigma && ("theta" %in% which)) {
    tr <- mmrm2sigma(x)
    sigmapar <- tr(coef(res))
    nam <- names(sigmapar)
    means <- if ("beta" %in% which) subset(res, seq_along(coef(x)))
    vars <- transform(res, tr)
    if (is.null(nam)) {
      vars <- labels(vars, paste0("sigma", seq_along(sigmapar)))
    }
    res <- vars
    if (!is.null(means)) {
      res <- c(means, vars)
    }
  }
  ## nam <- names(tr(coef(est)))
  res <- lava::estimate(res, ...)
  res$fit <- x
  res$sigma <- sigma
  structure(res, class=c("estimate.mmrm", "estimate"))
}

#' Conversion between covariance matrices and their upper-triangular elements
#'
#' `vec2sigma` constructs symmetric (covariance) matrices from a vector of
#' upper-triangular elements (column-major, including the diagonal), e.g., the
#' covariance parameters returned by `estimate(fit, sigma = TRUE)` for an mmrm
#' model. The vector may contain the elements of several matrices of the same
#' dimension (one per group), stacked after each other. `sigma2vec` is the
#' inverse operation.
#' @title Covariance matrix from upper-triangular elements
#' @param x (numeric) for `vec2sigma`, the vector of upper-triangular elements
#'   of \eqn{G} \eqn{p\times p} matrices (length \eqn{G p (p + 1) / 2}). For
#'   `sigma2vec`, a matrix or a (named) list of matrices.
#' @param groups (character) group names of the covariance matrices. If
#'   `NULL`, the groups are derived from the names of `x`, which must be of
#'   the form `<group><index>` (e.g., `grpA1, ..., grpAk, grpB1, ..., grpBk`).
#' @param visits (character) optional row and column names of the matrices.
#' @param simplify (logical) if `TRUE` and there is only a single group, the
#'   matrix is returned instead of a list.
#' @return `vec2sigma`: list of matrices named by group (or a matrix if
#'   `simplify = TRUE` and there is a single group). `sigma2vec`: numeric
#'   vector with the upper-triangular elements of each matrix.
#' @export
#' @examples
#' S <- matrix(c(2, 1, 1, 3), 2)
#' x <- sigma2vec(list(sigma = S))
#' x
#' vec2sigma(x)
#' vec2sigma(x, groups = "sigma", visits = c("v1", "v2"), simplify = TRUE)
#'
#' ## Two groups
#' vec2sigma(c(1, 0.5, 2, 3, 1, 4), groups = c("A", "B"))
vec2sigma <- function(x, groups=NULL, visits=NULL, simplify=FALSE) {
  if (is.null(groups)) { # derive groups from parameter names
    lbl <- names(x) # labeled as groupA1, ... groupAk, groupB1, ..., groupBK
    m <- regexpr("^.*?(?=[0-9]+$)", lbl, perl=TRUE) # match string###
    groups <- regmatches(lbl, m)
    ngroups <- length(unique(groups)) # groupA, groupB, ...
  } else {
    ngroups <- length(groups)
  }
  k <- length(x) / ngroups # k is the number of parameters in the upper-tri mat.
  p <- round((-1+sqrt(1+8*k))/2) # k = p*(p+1)/2 where sigma is a pxp matrix
  res <- c()
  for (i in seq_len(ngroups)) {
    sigma <- matrix(0, p, p)
    cur <- x[seq_len(k) + (i-1)*k]
    sigma[upper.tri(sigma, diag=TRUE)] <- cur
    for (i in seq_len(p-1)) { # make symmetric matrix
      for (j in seq(i+1, p)) {
        sigma[j, i] <- sigma[i, j]
      }
    }
    dimnames(sigma) <- list(visits, visits)
    res <- c(res, list(sigma))
  }
  names(res) <- unique(groups)
  if (simplify && length(res)==1L) return(res[[1]])
  return(res)
}

#' @rdname vec2sigma
#' @export
sigma2vec <- function(x) {
  if (is.matrix(x)) x <- list(x)
  unlist(lapply(x, function(v) v[upper.tri(v, diag=TRUE)]))
}

mmrm2sigma <- function(object) {
  function(p = pars(object), vec = TRUE) {
    np1 <- length(object$beta_est)
    np2 <- length(object$theta_est)
    if (identical(length(p), np1+np2)) {
      p <- p[-seq_len(np1)]
    }
    V <- .mmrm_varcor(object, p)
    if (!vec) return(V)
    sigma2vec(V)
  }
}

## Reconstruct the VarCorr-style covariance (matrix or list-per-group) for a
## given theta.  When theta is NULL (default) the result equals VarCorr(fit)
## exactly.
## Scope: discrete-visit covariances (us, cs, toep, ad, ar1)
.mmrm_varcor <- function(fit, theta = NULL) {
  stopifnot(inherits(fit, "mmrm"))
  if (is.null(theta)) theta <- mmrm::component(fit, "theta_est")

  cs <- mmrm::as.cov_struct(fit$formula_parts$formula)
  if (identical(mmrm::component(fit, "cov_type"), "sp_exp")) {
    stop("varcor_at(): spatial covariance (sp_exp) is not supported.")
  }

  rep <- fit$tmb_object$report(theta) # TMB API
  L   <- rep$covariance_lower_chol
  ntimes  <- mmrm::component(fit, "n_timepoints")
  ngroups   <- mmrm::component(fit, "n_groups")
  visit_names <- levels(fit$tmb_data$full_frame[[cs$visits]])

  cov_list <- lapply(seq_len(ngroups), function(g) {
    Lg <- L[seq((g - 1L) * ntimes + 1L, g * ntimes), , drop = FALSE]
    Sg <- tcrossprod(Lg)
    dimnames(Sg) <- list(visit_names, visit_names)
    Sg
  })
  if (ngroups==1L) {
    return(cov_list[[1]])
  }
  names(cov_list) <- levels(fit$tmb_data$subject_groups)
  cov_list
}

## Return an id ordered list per subject with:
##   id      : subject label (in order)
##   rows    : integer row indices into full_frame / x_matrix / y_vector
##             (ordered by IDs and visit levels)
##   visits  : character vector of observed visit-factor levels
##             (or numeric coordinates for spatial structures)
##   X, y, w : per-subject design, response, weights (ordeded by visits)
##   Sigma   : observed-visit covariance matrix (ordered by vists)
##   group   : character group label (or NA_character_ if no grouping)
##             (ordered by visits)
##
## This implementation focuses on discrete-time covariances:
## us, cs, toep, ad, ar1
.mmrm_subjects <- function(fit, theta=NULL) {
  stopifnot(inherits(fit, "mmrm"))
  cs        <- mmrm::as.cov_struct(fit$formula_parts$formula)
  subj_var  <- mmrm::component(fit, "subject_var")
  visit_var <- cs$visits
  group_var <- if (length(cs$group) == 0L) NA_character_ else cs$group

  ff <- mmrm::component(fit, "full_frame")
  X  <- mmrm::component(fit, "x_matrix")
  y  <- mmrm::component(fit, "y_vector")
  w  <- ff[["(weights)"]]
  if (is.null(w)) w <- rep(1, length(y))

  subj_full <- ff[[subj_var]]
  visit_full <- ff[[visit_var]]
  visit_char <- as.character(visit_full)

  if (!is.null(theta)) {
    Sig <- .mmrm_varcor(fit, theta)
  } else {
    Sig <- mmrm::VarCorr(fit)
  }
  is_grouped <- is.list(Sig) && !is.matrix(Sig)
  if (is_grouped) {
    group_full <- as.character(ff[[group_var]])
  } else {
    group_full <- rep(NA_character_, nrow(ff))
  }

  ## Split rows by subject preserving encounter order
  subj_char <- as.character(subj_full) # ordered list of IDs
  first_seen <- !duplicated(subj_char)
  subj_order <- subj_char[first_seen] # subject IDs in order
  rows_by <- split(seq_along(subj_char), subj_char)[subj_order]
  lapply(subj_order, function(sid) {
    rr <- rows_by[[sid]]
    vv <- visit_char[rr]
    gg <- group_full[rr[1]]
    Sfull <- if (is_grouped) Sig[[gg]] else Sig
    Si <- Sfull[vv, vv, drop = FALSE]
    vn <- match(vv, rownames(Sfull))
    list(id      = sid,
         rows    = rr,
         visits  = vv,
         visitn  = vn,
         group   = gg,
         X       = X[rr, , drop = FALSE],
         y       = y[rr],
         w       = w[rr],
         Sigma   = Si)
  })
}

## Weighted covariance of a subject. mmrm models the covariance of subject i as
##   Sigma_i = W_i^{-1/2} Sigma W_i^{-1/2}, W_i = diag(w_i) i.e., cov(y_ij,
##   y_ik) = Sigma_jk / sqrt(w_ij * w_ik).
.mmrm_wsigma <- function(S, w) {
  if (all(w == 1)) return(S)
  sw <- sqrt(w)
  return(S / tcrossprod(sw)) # S_jk / sqrt(w_j * w_k) = W^{-1/2} S W^{-1/2}
}

## Per-subject beta score: n_subj x p_beta.
.mmrm_score_beta <- function(fit, beta=NULL, theta=NULL) {
  subj <- .mmrm_subjects(fit, theta = theta)
  n <- length(subj)
  if (is.null(beta)) beta <- mmrm::component(fit, "beta_est")
  p <- length(beta)
  Sinv <- lapply(subj, function(s) solve(.mmrm_wsigma(s$Sigma, s$w)))
  r    <- lapply(subj, function(s) s$y - as.numeric(s$X %*% beta))
  Sir  <- Map(function(Si, ri) Si %*% ri, Sinv, r)
  Ub <- matrix(0, n, p,
               dimnames = list(NULL, names(beta)))
  for (i in seq_len(n)) {
    Ub[i, ] <- as.numeric(crossprod(subj[[i]]$X, Sir[[i]]))
  }
  Ub
}

## Derivative of the covariance matrix with respect to theta.
##
## mmrm uses TMB automatic differentiation internally but does not expose the
## raw dSigma/dtheta. Differentiate the reported covariance instead. The
## result has one element per theta parameter; each element is either a
## covariance matrix or a named list of covariance matrices, one per group.
dSigma_dtheta <- function(fit,
                          theta = NULL
                          ) {
  stopifnot(inherits(fit, "mmrm"))

  if (is.null(theta)) {
    theta <- pars(fit, "theta")
  }
  n_group <- mmrm::component(fit, "n_groups")
  cov_struct <- mmrm::as.cov_struct(fit$formula_parts$formula)
  visit_names <- levels(fit$tmb_data$full_frame[[cov_struct$visits]])
  group_names <- if (n_group > 1L) {
    levels(fit$tmb_data$subject_groups)
  } else {
    "sigma"
  }
  jacobian <- numDeriv::jacobian(
                          function(theta) {
                            sigma <- .mmrm_varcor(fit, theta)
                            sigma2vec(sigma)
                          },
                          theta,
                          method = lava::lava.options()$Dmethod,
                        )
  derivative <- lapply(seq_along(theta), function(k) jacobian[, k])
  derivative <- lapply(derivative,
                       function(x) {
                         vec2sigma(x,
                                   groups = group_names,
                                   visits = visit_names,
                                   simplify=TRUE)
                         })
  return(derivative)
}


## Restrict a covariance derivative to one subject's observed visits
## (including the subject's weights, see .mmrm_wsigma).
.subject_block <- function(dSigma, subject) {
  if (is.list(dSigma) && !is.matrix(dSigma)) {
    dSigma <- dSigma[[subject$group]]
  }
  dS <- dSigma[subject$visits, subject$visits, drop = FALSE]
  return(.mmrm_wsigma(dS, subject$w))
}

## Per-subject score for theta.
##
##   U_ik = -1/2 { tr(Sigma_i^-1 dSigma_i/dtheta_k)
##                 - r_i' Sigma_i^-1 dSigma_i/dtheta_k Sigma_i^-1 r_i }
.mmrm_score_theta <- function(fit,
                              theta = NULL,
                              beta = NULL
                              ) {
  subj <- .mmrm_subjects(fit, theta = theta)
  n <- length(subj)

  if (is.null(beta)) beta <- mmrm::component(fit, "beta_est")
  if (is.null(theta)) theta <- mmrm::component(fit, "theta_est")
  q <- length(theta)

  Sinv <- lapply(subj, function(s) solve(.mmrm_wsigma(s$Sigma, s$w)))
  r <- lapply(subj, function(s) s$y - as.numeric(s$X %*% beta))
  Sir <- Map(function(Si, ri) as.numeric(Si %*% ri), Sinv, r)
  dSig <- dSigma_dtheta(fit, theta = theta)

  Ut <- matrix(0, n, q, dimnames = list(NULL, paste0("theta", seq_len(q))))
  for (i in seq_len(n)) {
    Si <- Sinv[[i]]
    Siri <- Sir[[i]]
    for (k in seq_len(q)) {
      dS_ik <- .subject_block(dSig[[k]], subj[[i]])
      tr_term <- sum(Si * dS_ik)
      quad <- as.numeric(crossprod(Siri, dS_ik %*% Siri))
      Ut[i, k] <- -0.5 * (tr_term - quad)
    }
  }

  return(Ut)
}
