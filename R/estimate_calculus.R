# ---- Merge, subset ------------------------------------------------------

##' @export
merge.estimate <- function(x,y,...,
                           id,
                           paired=FALSE,
                           labels=NULL,
                           keep=NULL,
                           subset=NULL,
                           regex=FALSE,
                           drop.ic = FALSE,
                           ignore.case=FALSE) {
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
    trans <- unlist(lapply(
      objects, function(x) !is.null(x[["back.transform"]])
    ))
    if (any(trans)) warning("back-transformation ignored (`back.transform`)")
    if (length(nai <- names(objects)=="NA")>0)
    names(objects)[which(nai)] <- ""
    if (!missing(subset)) {
      if (regex) {
        # TODO
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
      ## if (!missing(id) && !is.null(id)) {
      ##   warning(
      ##     "Argument 'id' is only applicable to objects with an influence function. ",
      ##     "It will be used only for objects with IC; ",
      ##     "cross-terms involving objects without IC will be NA."
      ##   )
      ## }
      npar       <- unlist(lapply(objects, function(x) length(coef(x, messages=0))))
      col_ends   <- cumsum(npar)
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

    if (!missing(id) && is.null(id)) { ## Independence between datasets in x,y,...
        nn <- unlist(lapply(
          objects,
          function(x) nrow(x$IC)
        ))
        cnn <- c(0, cumsum(nn))
        id <- list()
        for (i in seq_along(nn)) {
          id <- c(id, list(seq(nn[i]) + cnn[i]))
        }
    }
    if (missing(id)) {
      if (paired) { ## One-to-one dependence between observations in x,y,...
        id <- lapply(objects, function(x) {
          seq_len(NROW(x$IC))
        })
        } else {
            id <- lapply(objects, function(x) x$id)
        }
    } else {
        nn <- unlist(lapply(objects,function(x) NROW(IC(x))))
        if (length(id)==1 && is.logical(id)) {
            if (id) {
                if (any(nn[1]!=nn)) stop("Expected objects of the same size: ", paste(nn,collapse=","))
                id0 <- seq(nn[1]); id <- c()
                for (i in seq(length(nn))) id <- c(id,list(id0))
            } else {
                id <- c()
                N <- cumsum(c(0,nn))
                for (i in seq(length(nn))) id <- c(id,list(seq(nn[i])+N[i]))
            }
        }
        if (length(id)!=length(objects)) stop("Same number of id-elements as model objects expected")
        idlen <- unlist(lapply(id,length))
        if (!identical(idlen,nn)) stop("Wrong lengths of 'id': ", paste(idlen,collapse=","), "; ", paste(nn,collapse=","))
    }
    ids <- ic_all <- c(); count <- 0
    for (z in objects) {
        count <- count+1
        clidx <- NULL
        id0 <- id[[count]]
        icz <- IC(z)
        if (is.null(id0)) {
            id0 <- rownames(icz)
            if (is.null(id0)) stop("Need id for object number ", count)
        }
        if (!missing(subset)) icz <- icz[,subset,drop=FALSE]
        if (!lava.options()$cluster.index) {
            ic0 <- matrix(unlist(by(icz,id0,colSums)),byrow=TRUE,ncol=ncol(icz))
            ids <- c(ids, list(sort(unique(id0))))
        } else {
            if (!requireNamespace("mets",quietly=TRUE)) stop("'mets' package required")
            clidx <- mets::cluster.index(id0,mat=icz,return.all=TRUE)
            ic0 <- clidx$X
            ids <- c(ids, list(id0[as.vector(clidx$firstclustid)+1]))
        }
        ic0 <- ic0*NROW(ic0)/length(id0)
        ic_all <- c(ic_all, list(ic0))
    }
    id <- unique(unlist(ids))
    ic0 <- matrix(NA, nrow=length(id),ncol=length(coefs))
    model.index <- c()
    colpos <- 0
    for (i in seq(length(objects))) {
        relpos <- seq_along(coef(objects[[i]], messages=0))
        if (!missing(subset)) relpos <- seq_along(subset)
        ic0[match(ids[[i]], id), relpos + colpos] <- ic_all[[i]]
        ## midx <- objects[[i]]$model.index
        ## if (!is.null(midx)) {
        ##   midx <- lapply(midx, function(x) {
        ##     intersect(x, relpos) + colpos
        ##   })
        ## } else {
          midx <- list(relpos + colpos)
        ## }
        model.index <- c(model.index, midx)
        colpos <- colpos+tail(relpos,1)
    }
    rownames(ic0) <- id
    ## Rescale each column according to I(obs)/pr(obs)
    for (i in seq(NCOL(ic0))) {
      pr <- mean(!is.na(ic0[,i]))
      ic0[,i] <- ic0[,i]/pr
    }
    ic0[is.na(ic0)] <- 0
    res <- estimate.default(
      coef = coefs, stack = FALSE, data = NULL,
      IC = ic0, id = id, keep = keep
      )
    if (is.null(keep)) {
      res$model.index <- model.index
    }
    return(res)
}

##' @export
"%++%.estimate" <- function(x, ...) {
  merge(x, ...)
}

##' @export
"c.estimate" <- function(...) {
  args <- list(...)
  n_args <- length(args)
  # Handle names robustly
  arg_names <- names(args)
  # If names are NULL, create a vector of empty strings
  if (is.null(arg_names)) {
    arg_names <- character(n_args)
  }
  is_estimate <- unlist(lapply(args, function(x)
    inherits(x, c("estimate"))
    ))
  merge_args <- c("drop.ic", "paired")
  not_merge_arg <- which(arg_names %ni% merge_args)
  if (!all(is_estimate[not_merge_arg])) { # fallback to default concatenation
    cl <- class(args[[1]])
    class(args[[1]]) <- "list"
    return(do.call(c, args))
  }
  lab <- arg_names[not_merge_arg]
  arg_names[not_merge_arg] <- ""
  names(args) <- arg_names
  res <- do.call(merge, args)
  newlabels <- names(coef(res))
  if (!is.null(lab)) {
    idx <- which(lab != "")
    newlabels[idx] <- lab[idx]
    return(labels(res, newlabels))
  }
  return(res)
}

##' @export
subset.estimate <- function(x, keep, ...) {
  estimate(x, keep = keep, ...)
}

##' @export
"[.estimate" <- function(x, i, ...) {
  subset(x, i, ...)
}

##' @export
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

# ---- Trigonometric Functions --------------------------------------------

##' @export
sin.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- sin(p)
    structure(y, grad = diag(cos(p), nrow = length(p)))
  }, ...)
}

##' @export
cos.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- cos(p)
    structure(y, grad = diag(-sin(p), nrow = length(p)))
  }, ...)
}

##' @export
tan.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- tan(p)
    structure(y, grad = diag(1 / cos(p)^2, nrow = length(p)))
  }, ...)
}

# ---- Inverse Trigonometric Functions ------------------------------------

##' @export
asin.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- asin(p)
    structure(y, grad = diag(1 / sqrt(1 - p^2), nrow = length(p)))
  }, ...)
}

##' @export
acos.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- acos(p)
    structure(y, grad = diag(-1 / sqrt(1 - p^2), nrow = length(p)))
  }, ...)
}

##' @export
atan.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- atan(p)
    structure(y, grad = diag(1 / (1 + p^2), nrow = length(p)))
  }, ...)
}

# ---- Hyperbolic Functions -----------------------------------------------

##' @export
sinh.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- sinh(p)
    structure(y, grad = diag(cosh(p), nrow = length(p)))
  }, ...)
}

##' @export
cosh.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- cosh(p)
    structure(y, grad = diag(sinh(p), nrow = length(p)))
  }, ...)
}

##' @export
tanh.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- tanh(p)
    structure(y, grad = diag(1 / cosh(p)^2, nrow = length(p)))
  }, ...)
}

# ---- Inverse Hyperbolic Functions ---------------------------------------

##' @export
asinh.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- asinh(p)
    structure(y, grad = diag(1 / sqrt(p^2 + 1), nrow = length(p)))
  }, ...)
}

##' @export
acosh.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- acosh(p)
    structure(y, grad = diag(1 / sqrt(p^2 - 1), nrow = length(p)))
  }, ...)
}

##' @export
atanh.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- atanh(p)
    structure(y, grad = diag(1 / (1 - p^2), nrow = length(p)))
  }, ...)
}

# ---- Other Common Functions ---------------------------------------------

##' @export
log1p.estimate <- function(x, ...) {
  # log(1 + p) — more numerically stable than log(p+1) for small p
  estimate(x, function(p) {
    y <- log1p(p)
    structure(y, grad = diag(1 / (1 + p), nrow = length(p)))
  }, ...)
}

##' @export
expm1.estimate <- function(x, ...) {
  # exp(p) - 1 — more numerically stable than exp(p)-1 for small p
  estimate(x, function(p) {
    y <- expm1(p)
    structure(y, grad = diag(exp(p), nrow = length(p)))
  }, ...)
}

##' @export
log.estimate <- function(x, base = exp(1), ...) {
  estimate(x, function(p) {
    y <- log(p, base = base)
    structure(y,
              grad = diag(1 / (p * log(base)), nrow=length(p)))
  }, ...)
}

##' @export
exp.estimate <- function(x, ...) {
  estimate(x, function(p) {
    y <- exp(p)
    structure(y, grad = diag(y, nrow=length(p)))
  }, ...)
}

##' @export
sqrt.estimate <- function(x, ...) {
  estimate(x^.5, ...)
}

##' @export
sum.estimate <- function(x, ...) {
  estimate(x,
           function(p)
             structure(sum(p),
                       grad = matrix(1, nrow = 1, ncol = length(p))),
           ...)
}

 ##' @export
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

##' @export
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

##' @export
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

##' @export
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

##' @export
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

##' @export
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

##' @export
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

# ---- == ---- ------------------------------------------------------------

##' @export
"==.estimate" <- function(e1, e2) {
  if (!(is.numeric(e1) || is.numeric(e2))) stop("numeric comperator needed")
  null <- if (is.numeric(e1)) e1 else e2
  e <- if (is.numeric(e1)) e2 else e1
  estimate(e, null=null)
}
