# ---- Merge, subset ------------------------------------------------------


##' @export
"[.estimate" <- function(x, i, ...) {
  subset(x, i, ...)
}

##' @export
"c.estimate" <- function(...) {
  args <- list(...)
  labels <- names(args)
  if (!is.null(labels)) {
    names(args) <- NULL
    args$labels <- labels
  }
  res <- do.call(merge, args)
  res$model.index <- NULL
  return(res)
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
"^.estimate" <- function(e1, e2, ...) {
  if (!is.numeric(e2) && length(e2) != 1L) stop("exponent should be a scalar")
  estimate(e1, function(p) {
    y <- p**e2
    structure(y, grad = diag(e2*p**(e2-1), nrow=length(p)))
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

##' @export
"+.estimate" <- function(e1, e2, ...) {
  operator_estimate(
    e1, e2,
    function(x, y, x_const=FALSE, y_const=FALSE) {
      nx <- length(x)
      ny <- length(y)
      n <- max(nx, ny)
      grad <-
        if (y_const) { # y is constant: df/dx only
          if (nx > 1L || ny == 1L) { # ny = nx or ny = 1
            diag(1, ncol=nx, nrow=nx)
          } else { # ny > 1, nx = 1
            matrix(1, nrow = n, ncol = 1L)
          }
        } else if (x_const) { # x is constant: df/dy only
          if (ny > 1L || nx == 1L) { # ny = nx or nx = 1
            diag(1, ncol=ny, nrow=ny)
          } else { # ny > 1, nx = 1
            matrix(1, nrow = n, ncol = 1L)
          }
        } else {
          # both estimates: df/d(c(x,y)) = [I | I]
          D <- matrix(0, nrow = n, ncol = nx + ny)
          if (nx == 1L) {
            D[, 1] <- 1
            D[, 2:ncol(D)] <- diag(1, ny, ny)
            D
          } else if (ny == 1L) {
            cbind(diag(1, nrow=n, ncol=n), 1)
          } else {
            cbind(diag(1, nrow = n, ncol=n), diag(1, nrow = n, ncol=n))
          }
        }
      structure(
        x + y,
        grad = grad
      )
      }, ...)
}

##' @export
"-.estimate" <- function(e1, e2, ...) {
  if (missing(e2)) return(-1*e1)
  operator_estimate(
    e1, e2,
    function(x, y, x_const=FALSE, y_const=FALSE) {
      nx <- length(x)
      ny <- length(y)
      n <- max(nx, ny)
      grad <-
        if (y_const) { # y is constant: df/dx only
          if (nx > 1L || ny == 1L) { # ny = nx or ny = 1
            diag(1, ncol=nx, nrow=nx)
          } else { # ny > 1, nx = 1
            matrix(1, nrow = n, ncol = 1L)
          }
        } else if (x_const) { # x is constant: df/dy only
          if (ny > 1L || nx == 1L) { # ny = nx or nx = 1
            diag(-1, ncol=ny, nrow=ny)
          } else { # ny > 1, nx = 1
            matrix(-1, nrow = n, ncol = 1L)
          }
        } else {
          # both estimates: df/d(c(x,y)) = [I | I]
          D <- matrix(0, nrow = n, ncol = nx + ny)
          if (nx == 1L) {
            D[, 1] <- 1
            D[, 2:ncol(D)] <- diag(-1, ny, ny)
            D
          } else if (ny == 1L) {
            cbind(diag(1, nrow=n, ncol=n), -1)
          } else {
            cbind(diag(1, nrow = n, ncol=n), diag(-1, nrow = n, ncol=n))
          }
        }
      structure(
        x - y,
        grad = grad
      )
      }, ...)
}

##' @export
"*.estimate" <- function(e1, e2, ...) {
  operator_estimate(
    e1, e2,
    function(x, y, x_const=FALSE, y_const=FALSE) {
      nx <- length(x)
      ny <- length(y)
      n <- max(nx, ny)
      grad <-
        if (y_const) { # y is constant: df/dx = y only
          if (nx > 1L || ny == 1L) { # ny = nx or ny = 1
            diag(y, ncol=nx, nrow=nx)
          } else { # ny > 1, nx = 1
            matrix(y, nrow = n, ncol = 1L)
          }
        } else if (x_const) { # x is constant: df/dy = x only
          if (ny > 1L || nx == 1L) { # ny = nx or nx = 1
            diag(x, ncol=ny, nrow=ny)
          } else { # ny > 1, nx = 1
            matrix(x, nrow = n, ncol = 1L)
          }
        } else {
          # both estimates: df/d(c(x,y)) = [y | x]
          D <- matrix(0, nrow = n, ncol = nx + ny)
          if (nx == 1L) {
            D[, 1] <- y
            D[, 2:ncol(D)] <- diag(x, ny, ny)
            D
          } else if (ny == 1L) {
            cbind(diag(y, nrow=n, ncol=n), x)
          } else {
            cbind(diag(y, nrow = n, ncol=n), diag(x, nrow = n, ncol=n))
          }
        }
      structure(
        x * y,
        grad = grad
      )
      }, ...)
}

##' @export
"/.estimate" <- function(e1, e2, ...) {
  operator_estimate(
    e1, e2,
    function(x, y, x_const=FALSE, y_const=FALSE) {
      nx <- length(x)
      ny <- length(y)
      n <- max(nx, ny)
      grad <-
        if (y_const) { # y is constant: df/dx = 1/y only
          if (nx > 1L || ny == 1L) { # ny = nx or ny = 1
            diag(1/y, ncol=nx, nrow=nx)
          } else { # ny > 1, nx = 1
            matrix(1/y, nrow = n, ncol = 1L)
          }
        } else if (x_const) { # x is constant: df/dy = -x/y^2 only
          if (ny > 1L || nx == 1L) { # ny = nx or nx = 1
            diag(-x / y^2, ncol=ny, nrow=ny)
          } else { # ny > 1, nx = 1
            matrix(-x / y^2, nrow = n, ncol = 1L)
          }
        } else {
          # both estimates: df/d(c(x,y)) = [1/y | -x/y^2]
          D <- matrix(0, nrow = n, ncol = nx + ny)
          if (nx == 1L) {
            D[, 1] <- 1/y
            D[, 2:ncol(D)] <- diag(-x / y^2, ny, ny)
            D
          } else if (ny == 1L) {
            cbind(diag(1/y, nrow=n, ncol=n), -x/y^2)
          } else {
            cbind(diag(1/y, nrow = n, ncol=n), diag(-x/y^2, nrow = n, ncol=n))
          }
        }
      structure(
        x / y,
        grad = grad
      )
      }, ...)
}
