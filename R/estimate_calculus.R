# ---- Merge, subset ------------------------------------------------------

#' Merge estimate objects
#'
#' @param x Object of class `estimate`
#' @param y Object of class `estimate`
#' @param ... Additional `estimate` objects or arguments
#' @param id Optional cluster variable
#' @param paired If TRUE a paired (matched) analysis is performed
#' @param labels Optional character vector of labels for the merged estimates
#' @param keep Optional character vector of parameter names to keep
#' @param subset Optional character vector of parameter names to subset
#' @param regex If TRUE, `keep` and `subset` are treated as regular expressions
#' @param sep Separator used for labeling
#' @param drop.ic If TRUE, drop the influence function from the result
#' @param ignore.case If TRUE, case is ignored in `keep`/`subset` matching
#' @param sort if true the returned influence function will be sorted according
#'   to the id variables (lexigraphically)
#' @seealso [c.estimate()]
#' @return Object of class `estimate` (see [estimate.default]).
#' @export
merge.estimate <- function(x, y,
                           ...,
                           id,
                           paired = FALSE,
                           labels = NULL,
                           keep = NULL,
                           subset = NULL,
                           regex = FALSE,
                           sep = FALSE,
                           drop.ic = FALSE,
                           ignore.case = FALSE,
                           sort = FALSE) {
    if (missing(y)) {
      objects <- c(list(x), list(...))
    } else {
      objects <- c(list(x), list(y), list(...))
    }
    if (drop.ic) {
      for (i in seq_along(objects))
      if (inherits(objects[[i]], "estimate")) {
        objects[[i]]$IC <- NULL
      }
    }
    if (length(nai <- names(objects)=="NA")>0)
    names(objects)[which(nai)] <- ""
    if (!missing(subset)) {
      if (regex) {
        warning("regular expression not supported by `merge` operation")
      }
      coefs <- unlist(lapply(objects, function(x) coef(x,messages=0)[subset]))
    } else {
      coefs <- unlist(lapply(objects, function(x) coef(x, messages=0)))
    }
    if (!is.null(labels)) {
      names(coefs) <- labels
    } else {
      names(coefs) <- make.unique(names(coefs))
    }
    if (regex) {
      if (!is.null(keep)) {
        cc <- names(coefs)
        keep <- unlist(lapply(keep, function(x) {
          cc[grepl(x, cc, perl = TRUE, ignore.case=ignore.case)]
        }))
      }
    }
    hasIC <- unlist(lapply(objects, function(x) !is.null(IC(x))))
    if (sum(!hasIC) > 0L) { # Some objects do not have influence function
      npar <- unlist(lapply(objects, function(x) length(coef(x, messages=0))))
      col_ends <- cumsum(npar)
      col_starts <- col_ends - npar + 1L
      V <- matrix(NA, nrow=sum(npar), ncol=sum(npar))

      ic_idx <- which(hasIC)
      if (length(ic_idx) > 0L) {
        if (length(ic_idx) > 1L) {
          ic_args <- c(objects[ic_idx],
                       list(paired = paired),
                       if (!missing(id)) list(id = id[ic_idx]))
          m_ic  <- do.call(merge, ic_args)
        } else m_ic <- objects[[ic_idx]]
        ic_cols <- unlist(Map(`:`, col_starts[ic_idx], col_ends[ic_idx]))
        V[ic_cols, ic_cols] <- vcov(m_ic)
      }
      # Fill diagonal blocks for non-IC objects
      for (k in which(!hasIC)) {
        pos <- col_starts[k]:col_ends[k]
        V[pos, pos] <- suppressMessages(vcov(objects[[k]]))
      }
      return(estimate(coef=coefs, vcov=V, keep=keep))
    }
    id_missing <- missing(id)
    id <- merge_ids(objects, id = if (!id_missing) id,
                    paired = paired, id_missing = id_missing)
    ics <- list(); model.index <- list(); colpos <- 0
    for (i in seq_along(objects)) {
        icz <- IC(objects[[i]])
        if (!missing(subset)) icz <- icz[, subset, drop = FALSE]
        if (lava.options()$check.ic) {
          check_ic_mean_zero(icz)
        }
        ics <- c(ics, list(cluster_sum_ic(icz, id[[i]])))
        model.index <- c(model.index, list(colpos + seq_len(ncol(icz))))
        colpos <- colpos + ncol(icz)
    }
    al <- align_ic(ics, lapply(id, unique))
    ic0 <- do.call(cbind, al$ic)
    id <- al$id
    if (sort) { # sort according to the ids (default: order of first IC)
      ord <- order(id)
      ic0 <- ic0[ord, , drop = FALSE]
      id <- id[ord]
    }
    rownames(ic0) <- id
    res <- estimate.default(
      coef = coefs, stack = FALSE, data = NULL,
      IC = ic0, id = id, keep = keep
      )
    if (is.null(keep) && sep) {
      res$model.index <- model.index
    }
    return(res)
}


## Ids of the rows of the influence functions of the estimate objects in
## 'objects' (used by merge.estimate):
##  - id=NULL or id=FALSE: independence (distinct ids across objects)
##  - id=TRUE or paired=TRUE: one-to-one matching (objects of the same size)
##  - 'id' missing (id_missing=TRUE): ids of the objects (index or rownames)
##  - otherwise a list of ids, one element for each object
merge_ids <- function(objects, id, paired = FALSE, id_missing = FALSE) {
  nn <- unlist(lapply(objects, function(x) NROW(IC(x))))
  if (!id_missing && (is.null(id) || isFALSE(id))) {
    cnn <- c(0, cumsum(nn))
    return(lapply(seq_along(nn), function(i) seq_len(nn[i]) + cnn[i]))
  }
  if ((id_missing && paired) || isTRUE(id)) {
    if (any(nn[1] != nn)) {
      stop("Expected objects of the same size: ", paste(nn, collapse = ","))
    }
    return(rep(list(seq_len(nn[1])), length(nn)))
  }
  if (id_missing) {
    return(lapply(seq_along(objects), function(i) {
      id0 <- ic_ids(objects[[i]], IC(objects[[i]]))
      if (is.null(id0)) stop("Need id for object number ", i)
      id0
    }))
  }
  if (length(id) != length(objects)) {
    stop("Same number of id-elements as model objects expected")
  }
  idlen <- unlist(lapply(id, length))
  if (!identical(idlen, nn)) {
    stop("Wrong lengths of 'id': ",
         paste(idlen, collapse = ","), "; ", paste(nn, collapse = ","))
  }
  id
}


#' @export
"%++%.estimate" <- function(x, ...) {
  merge(x, ...)
}

#' Concatenate estimate objects
#'
#' When all arguments are `estimate` objects, they are merged into a single
#' `estimate` object. When some arguments are not `estimate` objects but are
#' named numeric scalars/vectors, the merged `estimate` object is returned
#' with an `"extra"` attribute containing the provided values. This is useful
#' for bundling auxiliary per-iteration information (e.g., convergence status)
#' alongside an estimate for use with [sim.default()].
#'
#' @param ... `estimate` objects and/or named numeric values
#' @param as.list if TRUE the returned object will be of class `list` and not
#'   an `estimate` object.
#' @return An `estimate` object. If extra named numeric values were provided,
#'   the result carries an `"extra"` attribute (a named numeric vector).
#' @seealso [sim.default()] [merge.estimate()]
#' @details arguments `drop.ic`, `paired`, `sep` are passed to [merge.estimate]
#' @examples
#' e <- estimate(coef = c(a = 1, b = 2), vcov = diag(2) * 0.1)
#' # Bundle estimate with extra information
#' res <- c(e, converged = 1, niter = 10)
#' attr(res, "extra")
#' @export
c.estimate <- function(..., as.list = FALSE) {
  args <- list(...)
  if (as.list) { # fallback to default concatenation
    class(args[[1]]) <- "list"
    return(do.call(c, args))
  }
  # Handle names robustly
  arg_names <- names(args) %||% character(length(args))
  is_estimate <- vapply(args, inherits, logical(1L), "estimate")
  merge_args <- c("drop.ic", "paired", "sep")

  not_merge_arg <- which(arg_names %ni% merge_args)
  est_idx <- which(is_estimate)
  extra_idx <- setdiff(not_merge_arg, est_idx)
  lab <- arg_names[est_idx]
  # Collect extras from existing estimate objects (propagate "extra" attribute)
  extras_existing <- unlist(
    lapply(args[est_idx], function(x) attr(x, "extra"))
  )
  # Extract extra (non-estimate, non-merge) arguments
  extra <- c(extras_existing,
             if (length(extra_idx) > 0L) unlist(args[extra_idx]))
  # Blank non-merge names (estimate labels handled later, extras removed)
  arg_names[not_merge_arg] <- ""
  names(args) <- arg_names
  args[extra_idx] <- NULL
  # Merge estimate objects and apply potential merge_args
  res <- do.call(merge, args)
  # Add new labels
  newlabels <- names(coef(res))
  if (!is.null(lab)) {
    idx <- which(lab != "")
    newlabels[idx] <- lab[idx]
    res <- labels(res, newlabels)
  }
  # Attach extra as attribute (same pattern as c.summary.estimate)
  if (length(extra) > 0L) attr(res, "extra") <- extra
  res
}

#' @export
subset.estimate <- function(x, keep, ...) {
  estimate(x, keep = keep, ...)
}

#' @export
"[.estimate" <- function(x, i, ...) {
  subset(x, i, ...)
}

#' @export
with.estimate <- function(data, expr, ...) {
    # Recursively walk the expression tree and replace symbols
  # that match names in `data` with data["symbol"] calls
  replace_syms <- function(e) {
    # Base case: if it's a symbol, check if it matches a name in data
    if (is.symbol(e)) {
      nm <- as.character(e)
      if (nm %in% names(coef(data))) {
        # Replace symbol with data["nm"] call
        return(call("[", quote(data), nm))
      }
      return(e)
    }
    # Recursive case: walk the call tree
    if (is.call(e)) {
      return(as.call(lapply(e, replace_syms)))
    }
    # Literals (numbers, strings, etc.) — return as-is
    return(e)
  }
  # Substitute and transform the expression
  expr_sub  <- substitute(expr)
  expr_new  <- replace_syms(expr_sub)
  # Create a local environment where `data` exists,
  # with the parent frame as the enclosing environment
  eval_env <- new.env(parent = parent.frame())
  eval_env$data <- data
  # Evaluate in calling environment so non-estimate symbols still resolve
  eval(expr_new, envir = eval_env)
}

# ---- id / cluster  ------------------------------------------------------

#' @export
index.estimate <- function(x, ...) {
  return(x[["id"]])
}

#' @export
`index<-.estimate` <- function(x, ..., value) {
  estimate(x, id=value, ...)
}

# ---- Trigonometric Functions --------------------------------------------

#' @export
sin.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- sin(p)
    structure(y, grad = diag(cos(p), nrow = length(p)))
  }, ...)
}

#' @export
cos.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- cos(p)
    structure(y, grad = diag(-sin(p), nrow = length(p)))
  }, ...)
}

#' @export
tan.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- tan(p)
    structure(y, grad = diag(1 / cos(p)^2, nrow = length(p)))
  }, ...)
}

# ---- Inverse Trigonometric Functions ------------------------------------

#' @export
asin.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- asin(p)
    structure(y, grad = diag(1 / sqrt(1 - p^2), nrow = length(p)))
  }, ...)
}

#' @export
acos.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- acos(p)
    structure(y, grad = diag(-1 / sqrt(1 - p^2), nrow = length(p)))
  }, ...)
}

#' @export
atan.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- atan(p)
    structure(y, grad = diag(1 / (1 + p^2), nrow = length(p)))
  }, ...)
}

# ---- Hyperbolic Functions -----------------------------------------------

#' @export
sinh.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- sinh(p)
    structure(y, grad = diag(cosh(p), nrow = length(p)))
  }, ...)
}

#' @export
cosh.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- cosh(p)
    structure(y, grad = diag(sinh(p), nrow = length(p)))
  }, ...)
}

#' @export
tanh.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- tanh(p)
    structure(y, grad = diag(1 / cosh(p)^2, nrow = length(p)))
  }, ...)
}

# ---- Inverse Hyperbolic Functions ---------------------------------------

#' @export
asinh.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- asinh(p)
    structure(y, grad = diag(1 / sqrt(p^2 + 1), nrow = length(p)))
  }, ...)
}

#' @export
acosh.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- acosh(p)
    structure(y, grad = diag(1 / sqrt(p^2 - 1), nrow = length(p)))
  }, ...)
}

#' @export
atanh.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- atanh(p)
    structure(y, grad = diag(1 / (1 - p^2), nrow = length(p)))
  }, ...)
}

# ---- Other Common Functions ---------------------------------------------

#' @export
log1p.estimate <- function(x, ...) {
  # log(1 + p) — more numerically stable than log(p+1) for small p
  estimate(x, function(p) {
    y <- log1p(p)
    structure(y, grad = diag(1 / (1 + p), nrow = length(p)))
  }, ...)
}

#' @export
expm1.estimate <- function(x, ...) {
  # exp(p) - 1 — more numerically stable than exp(p)-1 for small p
  estimate(x, function(p) {
    y <- expm1(p)
    structure(y, grad = diag(exp(p), nrow = length(p)))
  }, ...)
}

#' @export
log.estimate <- function(x, base = exp(1), ...) {
  estimate(x, function(p) {
    y <- log(p, base = base)
    structure(y,
              grad = diag(1 / (p * log(base)), nrow=length(p)))
  }, ...)
}

#' @export
exp.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- exp(p)
    structure(y, grad = diag(y, nrow=length(p)))
  }, ...)
}

#' @export
sqrt.estimate <- function(x, ...) {
  estimate(x^.5, ...)
}

#' @export
sum.estimate <- function(x, ...) {
  estimate(x,
           function(p)
             structure(sum(p),
                       grad = matrix(1, nrow = 1, ncol = length(p))),
           ...)
}

 #' @export
"%*%.estimate" <- function(x, y, ...) {
  if (is.matrix(x)) {
    return(estimate(y, f=x, ...))
  } else if (is.matrix(y)) {
    return(estimate(x, f=t(y), ...))
  }
  sum(x * y)
}

prod_except <- function(p) {
  n     <- length(p)
  left  <- c(1, cumprod(p[-n]))
  right <- rev(cumprod(rev(p[-1])))
  right <- c(right, 1)
  left * right
}

#' @export
prod.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- prod(p)
    # grad of prod(p) is: d/dp_i = prod(p) / p_i
    structure(y, grad = matrix(prod_except(p), nrow = 1, ncol = length(p)))
  }, ...)
}

# ---- +,-,*,/ ------------------------------------------------------------

operator_estimate <- function(x, y, op, ...) {
  x_const <- is.numeric(x)
  y_const <- is.numeric(y)
  np1 <- ifelse(x_const, length(x), length(coef(x)))
  np2 <- ifelse(y_const, length(y), length(coef(y)))

  if (np1 != 1L && np2 != 1L && np1 != np2) {
    stop("expecting equal length objects or one of them to be a scalar")
  }
  if (y_const || x_const) { # x estimate, y numeric
    e <- if (y_const) x else y
    return(
      estimate(e, function(p) {
        if (y_const) {
          op(p, y, y_const=TRUE)
        } else {
          op(x, p, x_const=TRUE)
        }
      }, ...)
    )
  }
  e <- merge(x, y) # x estimate, y estimate
  estimate(e, function(p) {
    p1 <- p[seq_len(np1)]
    p2 <- p[seq_len(np2)+np1]
    res <- op(p1, p2)
    return(res)
  }, ...)
}

operator_grad <- function(x, y, x_const, y_const, dx, dy) {
  nx <- length(x)
  ny <- length(y)
  n <- max(nx, ny)
  grad <-
    if (y_const) { # y is constant: df/dx only
      if (nx > 1L || ny == 1L) { # ny = nx or ny = 1
        diag(dx, ncol=nx, nrow=nx)
      } else { # ny > 1, nx = 1
        matrix(dx, nrow = n, ncol = 1L)
      }
    } else if (x_const) { # x is constant: df/dy only
      if (ny > 1L || nx == 1L) { # ny = nx or nx = 1
        diag(dy, ncol=ny, nrow=ny)
      } else { # ny > 1, nx = 1
        matrix(dy, nrow = n, ncol = 1L)
      }
    } else {
      # both estimates: df/d(c(x,y)) = [dx | dy]
      D <- matrix(0, nrow = n, ncol = nx + ny)
      if (nx == 1L) {
        D[, 1] <- dx
        D[, 2:ncol(D)] <- diag(dy, ny, ny)
        D
      } else if (ny == 1L) {
        cbind(diag(dx, nrow=n, ncol=n), dy)
      } else {
        cbind(diag(dx, nrow = n, ncol=n), diag(dy, nrow = n, ncol=n))
      }
    }
  return(grad)
}

#' @export
"+.estimate" <- function(e1, e2, ...) {
  operator_estimate(
    e1, e2,
    function(x, y, x_const=FALSE, y_const=FALSE) {
      structure(
        x + y,
        grad = operator_grad(x, y, x_const, y_const,
                             dx = 1, dy = 1)
      )
      }, ...)
}

#' @export
"-.estimate" <- function(e1, e2, ...) {
  if (missing(e2)) return(-1*e1)
  operator_estimate(
    e1, e2,
    function(x, y, x_const=FALSE, y_const=FALSE) {
      structure(
        x - y,
        grad = operator_grad(x, y, x_const, y_const,
                             dx=1, dy=-1)
      )
      }, ...)
}

#' @export
"*.estimate" <- function(e1, e2, ...) {
  operator_estimate(
    e1, e2,
    function(x, y, x_const=FALSE, y_const=FALSE) {
      structure(
        x * y,
        grad = operator_grad(x, y, x_const, y_const,
                             dx=y, dy=x)
      )
      }, ...)
}

#' @export
"/.estimate" <- function(e1, e2, ...) {
  operator_estimate(
    e1, e2,
    function(x, y, x_const=FALSE, y_const=FALSE) {
      structure(
        x / y,
        grad = operator_grad(x, y, x_const, y_const,
                             dx = 1 / y, dy = -x / y^2)
      )
      }, ...)
}

#' @export
"^.estimate" <- function(e1, e2, ...) {
  operator_estimate(
    e1, e2,
    function(x, y, x_const=FALSE, y_const=FALSE) {
      nx <- length(x)
      ny <- length(y)
      n <- max(nx, ny)
      f <- x^y
      structure(
        f,
        grad = operator_grad(x, y, x_const, y_const,
                             dx = y*x^(y-1), dy = log(x)*f)
      )
    }, ...)
}

# ---- == / hypothesis ----------------------------------------------------

#' @export
"==.estimate" <- function(e1, e2) {
  if (!(is.numeric(e1) || is.numeric(e2))) stop("numeric comparator needed")
  null <- if (is.numeric(e1)) e1 else e2
  e <- if (is.numeric(e1)) e2 else e1
  summary(e, null = null)
}
