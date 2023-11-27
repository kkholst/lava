
##' @export
estimate <- function(x, ...) UseMethod("estimate")

##' @export
estimate.list <- function(x, ...) {
  if (inherits(x[[1]], "lvm")) return(estimate.lvmlist(x, ...))
  lapply(x, function(x) estimate(x, ...))
}

##' @export
estimate.data.frame <- function(x, ...) {
  estimate(as.matrix(x), ...)
}


IC.quantile <- function(x, estimate, probs=0.5, ...) {
  x <- na.omit(x)
  f0 <- density(x, ...)
  ## U <- function(est) (tau - (x <= est))
  if (missing(estimate)) {
    estimate <- quantile(x, probs=probs)
  }
  res <- c()
  for (i in seq_len(length(estimate))) {
    res <- cbind(res,
    (probs[i] - (x <= estimate[i]))/with(f0, approx(x, y, estimate[i]))$y
    )
  }
  res
}

##' @export
estimate.array <- function(x, type="mean", probs=0.5, ...) {
  if (missing(x) || is.null(x)) {
    return(estimate(NULL, ...))
  }
  dots <- list(...)
  density.args <- dots[]
  cc <- apply(x, 2, function(y) mean(y, na.rm = TRUE))
  ic <- apply(x, 2, function(y) y - mean(y, na.rm = TRUE))
  if (tolower(type) %in% c("var", "variance")) {
    n <- NROW(x)
     cc <- apply(x, 2, function(y) mean((y - mean(y)^2), na.rm = TRUE))
    ic <- ic^2
    for (i in seq_len(NCOL(ic))) {
      ic[, i] <- ic[, i] - cc[i]
    }
  }
  if (tolower(type) %in% c("quantile")) {
    density.args <- list()
    dargs <- names(formals(density.default))
    didx <- which(dargs %in% names(dots))
    if (length(didx)>0) {
      density.args <- dots[dargs[didx]]
      dots[dargs[didx]] <- NULL
    }
    cc <- unlist(apply(x, 2, function(y)
      quantile(y, probs=probs, na.rm = TRUE),
      simplify=FALSE))
    ic <- c()
    for (i in seq_len(NCOL(x))) {
      ic <- cbind(ic, do.call(IC.quantile,
                        c(list(x[,i], probs=probs), density.args)))
    }
  }
  if (any(c("vcov", "IC") %in% names(list(...)))) {
    return(estimate(NULL, coef = cc, ...))
  }
  do.call(estimate, c(list(NULL, coef=cc, IC=ic), dots))
}

##' Estimation of functional of parameters
##'
##' Estimation of functional of parameters.
##' Wald tests, robust standard errors, cluster robust standard errors,
##' LRT (when \code{f} is not a function)...
##' @param x model object (\code{glm}, \code{lvmfit}, ...)
##' @param f transformation of model parameters and (optionally) data, or contrast matrix (or vector)
##' @param ... additional arguments to lower level functions
##' @param data \code{data.frame}
##' @param id (optional) id-variable corresponding to ic decomposition of model parameters.
##' @param iddata (optional) id-variable for 'data'
##' @param stack if TRUE (default)  the i.i.d. decomposition is automatically stacked according to 'id'
##' @param average if TRUE averages are calculated
##' @param subset (optional) subset of data.frame on which to condition (logical expression or variable name)
##' @param score.deriv (optional) derivative of mean score function
##' @param level level of confidence limits
##' @param IC if TRUE (default) the influence function decompositions are also returned (extract with \code{IC} method)
##' @param type type of small-sample correction
##' @param keep (optional) index of parameters to keep from final result
##' @param use (optional) index of parameters to use in calculations
##' @param regex If TRUE use regular expression (perl compatible) for keep,use arguments
##' @param ignore.case Ignore case-sensitiveness in regular expression
##' @param contrast (optional) Contrast matrix for final Wald test
##' @param null (optional) null hypothesis to test
##' @param vcov (optional) covariance matrix of parameter estimates (e.g. Wald-test)
##' @param coef (optional) parameter coefficient
##' @param robust if TRUE robust standard errors are calculated. If
##' FALSE p-values for linear models are calculated from t-distribution
##' @param df degrees of freedom (default obtained from 'df.residual')
##' @param print (optional) print function
##' @param labels (optional) names of coefficients
##' @param label.width (optional) max width of labels
##' @param only.coef if TRUE only the coefficient matrix is return
##' @param back.transform (optional) transform of parameters and confidence intervals
##' @param folds (optional) aggregate influence functions (divide and conquer)
##' @param cluster (obsolete) alias for 'id'.
##' @param R Number of simulations (simulated p-values)
##' @param null.sim Mean under the null for simulations
##' @details
##'
##' influence function decomposition of estimator \eqn{\widehat{\theta}} based on
##' data \eqn{Z_1,\ldots,Z_n}:
##' \deqn{\sqrt{n}(\widehat{\theta}-\theta) = \frac{1}{\sqrt{n}}\sum_{i=1}^n IC(Z_i; P) + o_p(1)}
##' can be extracted with the \code{IC} method.
##'
##' @export
##' @export estimate.default
##' @examples
##'
##' ## Simulation from logistic regression model
##' m <- lvm(y~x+z);
##' distribution(m,y~x) <- binomial.lvm("logit")
##' d <- sim(m,1000)
##' g <- glm(y~z+x,data=d,family=binomial())
##' g0 <- glm(y~1,data=d,family=binomial())
##'
##' ## LRT
##' estimate(g,g0)
##'
##' ## Plain estimates (robust standard errors)
##' estimate(g)
##'
##' ## Testing contrasts
##' estimate(g,null=0)
##' estimate(g,rbind(c(1,1,0),c(1,0,2)))
##' estimate(g,rbind(c(1,1,0),c(1,0,2)),null=c(1,2))
##' estimate(g,2:3) ## same as cbind(0,1,-1)
##' estimate(g,as.list(2:3)) ## same as rbind(c(0,1,0),c(0,0,1))
##' ## Alternative syntax
##' estimate(g,"z","z"-"x",2*"z"-3*"x")
##' estimate(g,z,z-x,2*z-3*x)
##' estimate(g,"?")  ## Wildcards
##' estimate(g,"*Int*","z")
##' estimate(g,"1","2"-"3",null=c(0,1))
##' estimate(g,2,3)
##'
##' ## Usual (non-robust) confidence intervals
##' estimate(g,robust=FALSE)
##'
##' ## Transformations
##' estimate(g,function(p) p[1]+p[2])
##'
##' ## Multiple parameters
##' e <- estimate(g,function(p) c(p[1]+p[2],p[1]*p[2]))
##' e
##' vcov(e)
##'
##' ## Label new parameters
##' estimate(g,function(p) list("a1"=p[1]+p[2],"b1"=p[1]*p[2]))
##' ##'
##' ## Multiple group
##' m <- lvm(y~x)
##' m <- baptize(m)
##' d2 <- d1 <- sim(m,50,seed=1)
##' e <- estimate(list(m,m),list(d1,d2))
##' estimate(e) ## Wrong
##' ee <- estimate(e, id=rep(seq(nrow(d1)), 2)) ## Clustered
##' ee
##' estimate(lm(y~x,d1))
##'
##' ## Marginalize
##' f <- function(p,data)
##'   list(p0=lava:::expit(p["(Intercept)"] + p["z"]*data[,"z"]),
##'        p1=lava:::expit(p["(Intercept)"] + p["x"] + p["z"]*data[,"z"]))
##' e <- estimate(g, f, average=TRUE)
##' e
##' estimate(e,diff)
##' estimate(e,cbind(1,1))
##'
##' ## Clusters and subset (conditional marginal effects)
##' d$id <- rep(seq(nrow(d)/4),each=4)
##' estimate(g,function(p,data)
##'          list(p0=lava:::expit(p[1] + p["z"]*data[,"z"])),
##'          subset=d$z>0, id=d$id, average=TRUE)
##'
##' ## More examples with clusters:
##' m <- lvm(c(y1,y2,y3)~u+x)
##' d <- sim(m,10)
##' l1 <- glm(y1~x,data=d)
##' l2 <- glm(y2~x,data=d)
##' l3 <- glm(y3~x,data=d)
##'
##' ## Some random id-numbers
##' id1 <- c(1,1,4,1,3,1,2,3,4,5)
##' id2 <- c(1,2,3,4,5,6,7,8,1,1)
##' id3 <- seq(10)
##'
##' ## Un-stacked and stacked i.i.d. decomposition
##' IC(estimate(l1,id=id1,stack=FALSE))
##' IC(estimate(l1,id=id1))
##'
##' ## Combined i.i.d. decomposition
##' e1 <- estimate(l1,id=id1)
##' e2 <- estimate(l2,id=id2)
##' e3 <- estimate(l3,id=id3)
##' (a2 <- merge(e1,e2,e3))
##'
##' ## If all models were estimated on the same data we could use the
##' ## syntax:
##' ## Reduce(merge,estimate(list(l1,l2,l3)))
##'
##' ## Same:
##' IC(a1 <- merge(l1,l2,l3,id=list(id1,id2,id3)))
##'
##' IC(merge(l1,l2,l3,id=TRUE)) # one-to-one (same clusters)
##' IC(merge(l1,l2,l3,id=FALSE)) # independence
##'
##'
##' ## Monte Carlo approach, simple trend test example
##'
##' m <- categorical(lvm(),~x,K=5)
##' regression(m,additive=TRUE) <- y~x
##' d <- simulate(m,100,seed=1,'y~x'=0.1)
##' l <- lm(y~-1+factor(x),data=d)
##'
##' f <- function(x) coef(lm(x~seq_along(x)))[2]
##' null <- rep(mean(coef(l)),length(coef(l))) ## just need to make sure we simulate under H0: slope=0
##' estimate(l,f,R=1e2,null.sim=null)
##'
##' estimate(l,f)
##' @aliases estimate estimate.default estimate.estimate merge.estimate estimate.array estimate.mlm
##' @method estimate default
##' @export
estimate.default <- function(x=NULL, f=NULL, ..., data, id,
                             iddata, stack=TRUE, average=FALSE, subset,
                             score.deriv, level=0.95, IC=robust,
                             type=c("robust", "df", "mbn"),
                             keep, use,
                             regex=FALSE, ignore.case=FALSE,
                             contrast, null, vcov, coef,
                             robust=TRUE, df=NULL,
                             print=NULL, labels, label.width,
                             only.coef=FALSE, back.transform=NULL,
                             folds=0,
                             cluster,
                             R=0,
                             null.sim) {
  cl <- match.call(expand.dots=TRUE)
  if ("iid" %in% names(cl)) {
    stop("The 'iid' argument is obsolete. Please use the 'IC' argument")
  }
  if (!missing(use)) {

    p0 <- c("f","contrast","only.coef","subset","average","keep","labels","null")
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
    pp <- suppressWarnings(try(stats::coef(x), "try-error"))
    if (inherits(x, "survreg") && length(pp) < NROW(x$var)) {
      pp <- c(pp, scale=x$scale)
    }
  }
  if (!missing(cluster)) id <- cluster
  if (expr || is.character(f) || (is.numeric(f) && !is.matrix(f))) { ## || is.call(f)) {
    dots <- lapply(substitute(placeholder(...))[-1], function(x) x)
    args <- c(list(
      coef = names(pp),
      x = substitute(f),
      regex = regex
    ), dots)
    f <- do.call(parsedesign, args)
  }
  if (!is.null(f) && !is.function(f)) {
    if (!(is.matrix(f) | is.vector(f)))
      return(compare(x, f, ...))
    contrast <- f
    f <- NULL
  }

  if (lava.options()$cluster.index) {
    if (!requireNamespace("mets", quietly=TRUE)) stop("'mets' package required")
  }
  if (missing(data))
    data <- tryCatch(model.frame(x), error=function(...) NULL)
  alpha <- 1 - level
  alpha.str <- paste(c(alpha/2, 1 -alpha/2)*100, "", sep="%")
  nn <- NULL
  if ((( (is.logical(IC) && IC) || length(IC)>0) && robust) &&
      (missing(vcov) || is.null(vcov) ||
       (is.logical(vcov) && vcov[1]==FALSE && !is.na(vcov[1])))) {
    ## If user supplied vcov, then don't estimate IC
    if (missing(score.deriv)) {
      if (!is.logical(IC)) {
        ic_theta <- cbind(IC)
        IC <- TRUE
      } else {
        suppressWarnings(ic_theta <- IC(x, folds=folds))
      }
    } else {
      suppressWarnings(ic_theta <- IC(x, score.deriv=score.deriv, folds=folds))
    }
  } else {
    if (!is.null(x) && (missing(vcov) ||
                        (is.logical(vcov) && !is.na(vcov)[1])))
      suppressWarnings(vcov <- stats::vcov(x))
    ic_theta <- NULL
  }

  if (any(is.na(ic_theta))) {
    ## Rescale each column according to I(obs)/pr(obs)
    for (i in seq(NCOL(ic_theta))) {
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
  if (missing(id) && inherits(x,"estimate") && !is.null(x$id))
    id <- x$id
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
      id <- if(is.null(ic_theta)) seq(nrow(data)) else seq(nprev)
      stack <- FALSE
    }
    if (is.character(id) && length(id)==1)
      id <- data[,id,drop=TRUE]
    if (!is.null(ic_theta)) {
      if (length(id)!=nprev) {
        if (!is.null(x$na.action) && (length(id)==length(x$na.action) + nprev)) {
          warning("Applying na.action")
          id <- id[-x$na.action]
        } else stop("Dimensions of i.i.d decomposition and 'id' does not agree")
      }
    } else {
      if (length(id)!=nrow(data)) {
        if (!is.null(x$na.action) && (length(id)==length(x$na.action)+nrow(data))) {
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
        ic_theta <- matrix(unlist(by(ic_theta,id,colSums)),
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
  if (!robust) {
    if (inherits(x,"lm") && family(x)$family=="gaussian"
        && is.null(df))
      df <- x$df.residual
    if (missing(vcov) && !is.null(x))
      suppressWarnings(vcov <- stats::vcov(x))
  }

  if (!is.null(ic_theta) && robust && (missing(vcov) || is.null(vcov))) {
    V <- var_ic(ic_theta)
    ## Small-sample corrections for clustered data
    K <- NROW(ic_theta)
    N <- attributes(ic_theta)$N
    if (is.null(N)) N <- K
    p <- NCOL(ic_theta)
    adj0 <- K/(K-p) ## Mancl & DeRouen, 2001
    adj1 <- K/(K-1) ## Mancl & DeRouen, 2001
    adj2 <- (N-1)/(N-p)*(K/(K-1)) ## Morel,Bokossa & Neerchal, 2003
    if (tolower(type[1])=="mbn" && !is.null(attributes(ic_theta)$bread)) {
      V0 <- V
      iI0 <- attributes(ic_theta)$bread
      I0 <- Inverse(iI0)
      delta <- min(0.5 ,p/(K-p))
      phi <- max(1, tr(I0%*%V0)*adj2/p)
      V <- adj2*V0 + delta*phi*iI0
    }
    if (tolower(type[1])=="df") {
      V <- adj0*V
    }
    if (tolower(type[1])=="df1") {
      V <- adj1*V
    }
    if (tolower(type[1])=="df2") {
      V <- adj2*V
    }
  } else {
    if (!missing(vcov)) {
      if (length(vcov)==1 && is.na(vcov)) vcov <- matrix(NA,length(pp),length(pp))
      V <- cbind(vcov)
    } else {
      suppressWarnings(V <- stats::vcov(x))
    }
  }

  ## Simulate p-value
  if (R>0) {
    if (is.null(f)) stop("Supply function 'f'")
    if (missing(null.sim)) null.sim <- rep(0, length(pp))
    est <- f(pp)
    if (is.list(est)) {
      nn <- names(est)
      est <- unlist(est)
      names(est) <- nn
    }
    if (missing(labels)) {
      labels <- colnames(rbind(est))
    }
    res <- simnull(R,f,mu=null.sim,sigma=V,labels=labels)
    return(structure(res, class=c("estimate.sim","sim"),
                     coef=pp,
                     vcov=V,
                     f=f,
                     estimate=est))
  }

  if (!is.null(f)) {
    form <- names(formals(f))
    dots <- ("..."%in%names(form))
    form0 <- setdiff(form, "...")
    parname <- "p"

    if (!is.null(form)) parname <- form[1] # unless .Primitive
    if (length(form0)==1 && !(form0%in%c("object", "data"))) {
      ##names(formals(f))[1] <- "p"
      parname <- form0
    }
    if (!is.null(ic_theta) ) {
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
      if (!average || (N<NROW(data))) { ## || NROW(data)==0)) { ## transformation not depending on data
        pp <- structure(as.vector(val), names=names(val))
        ic_theta <- ic_theta%*%t(D)
        V <- var_ic(ic_theta)
      } else {
        if (k>1) { ## More than one parameter (and depends on data)
          if (!missing(subset)) { ## Conditional estimate
            val <- apply(val,2,function(x) x*subset)
          }
          D0 <- matrix(nrow=k, ncol=length(pp))
          for (i in seq_len(k)) {
            D1 <- D[seq(N)+(i-1)*N, , drop=FALSE]
            if (!missing(subset)) ## Conditional estimate
              D1 <- apply(D1, 2, function(x) x*subset)
            D0[i,] <- colMeans(D1)
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
        ic1 <- (cbind(val)-rbind(pp)%x%cbind(rep(1,N)))
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
            message("Assuming independence between model iid decomposition and new data frame")
            V <- var_ic(ic1) + var_ic(ic2)
          } else {
            ic_theta <- ic1+ic2
            V <- var_ic(ic_theta)
          }
        }
      }
    }
  }

  if (is.null(V)) {
    res <- cbind(pp, NA, NA, NA, NA)
  } else {
    if (length(pp)==1)
      res <- rbind(c(pp, diag(V)^0.5))
    else
      res <- cbind(pp, diag(V)^0.5)
    beta0 <- res[, 1]

    if (!missing(null) && missing(contrast))
      beta0 <- beta0-null
    if (!is.null(df)) {
      za <- qt(1-alpha/2, df=df)
      pval <- 2*pt(abs(res[,1]/res[,2]), df=df, lower.tail=FALSE)
    } else {
      za <- qnorm(1-alpha/2)
      pval <- 2*pnorm(abs(res[,1]/res[,2]), lower.tail=FALSE)
    }
    res <- cbind(res, res[,1]-za*res[,2], res[,1]+za*res[,2], pval)
  }
  colnames(res) <- c("Estimate","Std.Err",alpha.str,"P-value")

  if (nrow(res)>0)
    if (!is.null(nn)) {
      rownames(res) <- nn
    } else {
      nn <- attributes(res)$varnames
      if (!is.null(nn))
        rownames(res) <- nn
      if (is.null(rownames(res)))
        rownames(res) <- paste0("p", seq(nrow(res)))
    }

  if (NROW(res)==0L) {
    coefs <- NULL
  } else {
    coefs <- res[, 1, drop=TRUE]
    names(coefs) <- rownames(res)
  }
  res <- structure(list(coef=coefs, coefmat=res, vcov=V,
                        IC=NULL, print=print, id=idstack),
                   class="estimate")
  if (IC) ## && is.null(back.transform))
    res$IC <- ic_theta
  if (length(coefs)==0L) return(res)

  if (!missing(contrast) | !missing(null)) {
    p <- length(res$coef)
    if (missing(contrast)) contrast <- diag(nrow=p)
    if (missing(null)) null <- 0
    if (is.vector(contrast) || is.list(contrast)) {
      contrast <- contr(contrast, names(res$coef), ...)
      ## if (length(contrast)==p) contrast <- rbind(contrast)
      ## else {
      ##     cont <- contrast
      ##     contrast <- diag(nrow=p)[cont,,drop=FALSE]
      ## }
    }
    cc <- compare(res, contrast=contrast, null=null,
                  vcov=V, level=level, df=df)
    res <- structure(c(res, list(compare=cc)), class="estimate")
    if (!is.null(df)) {
      pval <- with(cc, pt(abs(estimate[,1]-null)/estimate[, 2],
                          df=df, lower.tail=FALSE)*2)
    } else {
      pval <- with(cc, pnorm(abs(estimate[, 1]-null)/estimate[, 2],
                             lower.tail=FALSE)*2)
    }
    res$coefmat <- with(cc, cbind(estimate, pval))
    colnames(res$coefmat)[5] <- "P-value"
    rownames(res$coefmat) <- cc$cnames
    if (!is.null(res$IC)) {
      res$IC <- res$IC%*%t(contrast)
      colnames(res$IC) <- cc$cnames
    }
    res$compare$estimate <- NULL
    res$coef <- res$compare$coef
    res$vcov <- res$compare$vcov
    names(res$coef) <- gsub("(^\\[)|(\\]$)", "",
                            rownames(res$coefmat))
  }

  if (!is.null(back.transform)) {
    res$coefmat[, c(1, 3, 4)] <- do.call(back.transform,
                                         list(res$coefmat[,c(1, 3, 4)]))
    res$coefmat[, 2] <- NA
  }

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
  if (!missing(labels)) {
    names(res$coef) <- labels
    if (!is.null(res$IC))
      colnames(res$IC) <- labels
    if (!is.null(res$vcov))
      colnames(res$vcov) <- rownames(res$vcov) <- labels
    rownames(res$coefmat) <- labels
  }
  if (!missing(label.width)) {
    rownames(res$coefmat) <- make.unique(
      unlist(lapply(rownames(res$coefmat),
                    function(x) toString(x, width=label.width)))
    )
  }
  if (only.coef) return(res$coefmat)
  res$call <- cl
  res$back.transform <- back.transform
  res$n <- nrow(data)
  res$ncluster <- nrow(res$IC)
  return(res)
}

simnull <- function(R, f, mu, sigma, labels=NULL) {
  X <- rmvn0(R, mu=mu, sigma=sigma)
  est <- f(mu)
  res <- apply(X,1,f)
  if (is.list(est)) {
    nn <- names(est)
    est <- unlist(est)
    names(est) <- nn
    res <- matrix(unlist(res), byrow=TRUE, ncol=length(est))
  } else {
    res <- t(rbind(res))
  }
  if (is.null(labels)) {
    labels <- colnames(rbind(est))
    if (is.null(labels))
      labels <- paste0("p", seq_along(est))
  }
  colnames(res) <- labels
  return(res)
}

##' @export
estimate.estimate.sim <- function(x, f, R=0, labels, ...) {
  atr <- attributes(x)
  if (R>0) {
    if (missing(f)) {
      val <- simnull(R,f=atr[["f"]],mu=atr[["coef"]],sigma=atr[["vcov"]])
      res <- rbind(x,val)
      for (a in setdiff(names(atr),c("dim","dimnames")))
        attr(res,a) <- atr[[a]]
    } else {
      res <- simnull(R,f=f,mu=atr[["coef"]],sigma=atr[["vcov"]])
      for (a in setdiff(names(atr),c("dim","dimnames","f")))
        attr(res,a) <- atr[[a]]
      attr(f,"f") <- f
      est <- unlist(f(atr[["coef"]]))
      if (missing(labels)) labels <- colnames(rbind(est))
      attr(res,"estimate") <- est
    }
    if (!missing(labels)) colnames(res) <- labels
    return(res)
  }
  if (missing(f)) {
    if (!missing(labels)) colnames(res) <- labels
    return(x)
  }

  est <- f(atr[["coef"]])
  res <- apply(x,1,f)
  if (is.list(est)) {
    res <- matrix(unlist(res),byrow=TRUE,ncol=length(est))
  } else {
    res <- t(rbind(res))
  }
  if (missing(labels)) {
    labels <- colnames(rbind(est))
    if (is.null(labels)) labels <- paste0("p",seq_along(est))
  }
  colnames(res) <- labels
  for (a in setdiff(names(atr),c("dim","dimnames","f","estimate")))
    attr(res,a) <- atr[[a]]
  attr(f,"f") <- f
  attr(res,"estimate") <- unlist(est)
  return(res)
}

##' @export
print.estimate.sim <- function(x, level=.95, ...) {
  quantiles <- c((1-level)/2, 1-(1-level)/2)
  est <- attr(x, "estimate")
  mysummary <- function(x, INDEX, ...) {
    x <- as.vector(x)
    res <- c(mean(x, na.rm=TRUE),
             sd(x, na.rm=TRUE),
             quantile(x, quantiles, na.rm=TRUE),
             est[INDEX],
             mean(abs(x)>abs(est[INDEX]), na.rm=TRUE))

    names(res) <- c("Mean", "SD", paste0(quantiles*100, "%"),
                    "Estimate", "P-value")
    res
  }
  env <- new.env()
  assign("est", attr(x,"estimate"), env)
  environment(mysummary) <- env
  print(summary(x, fun=mysummary, ...))
}

##' @export
estimate.glm <- function(x, ...) {
  estimate.default(x, ...)
}

##' @export
IC.mlm <- function(x, ...) {
  cc <- coef(x)
  r <- residuals(x)
  X <- model.matrix(x)
  w <- weights(x)
  if (!is.null(w)) {
    r <- apply(r, 2, function(x) x * w)
  }
  p <- NROW(cc)
  q <- NCOL(cc)
  ics <- lapply(1:q, function(i) {
    ic <- apply(X, 2, function(x) x * r[, i])
  })
  res <- Reduce("cbind", ics)
  colnames(res) <- names(pars(x))
  return(res)
}

##' @export
pars.mlm <- function(x, ...) {
  cc <- coef(x)
  q <- NCOL(cc)
  nn <- unlist(lapply(
    1:q, function(i)
      paste0(colnames(cc)[i], ":", rownames(cc))
  ))
  coefs <- unlist(lapply(1:q, function(x) cc[, x, drop=TRUE]))
  names(coefs) <- nn
  coefs
}

##' @export
estimate.mlm <- function(x, ...) {
  estimate.default(x, coef=pars(x), ...)
}

##' @export
print.estimate <- function(x, type=0L, digits=4L, width=25L,
                           std.error=TRUE, p.value=TRUE,
                           sep=cli::symbol[["line"]], sep.which, sep.labels=NULL,
                           indent=" ", unique.names=TRUE,
                           na.print="", ...) {

  if (!is.null(x$print)) {
    x$print(x,digits=digits,width=width,...)
    return(invisible(x))
  }
  if (type>0 && !is.null(x$call)) {
    cat("Call: "); print(x$call)
    print(cli::rule())
  }
  if (type>0) {
    if (!is.null(x[["n"]]) && !is.null(x[["k"]])) {
      cat("n = ", x[["n"]], ", clusters = ", x[["k"]], "\n\n", sep="")
    } else {
      if (!is.null(x[["n"]])) {
        cat("n = ", x[["n"]], "\n\n", sep="")
      }
      if (!is.null(x[["k"]])) {
        cat("n = ", x[["k"]], "\n\n",sep="")
      }
    }
  }

  cc <- x$coefmat
  if (!is.null(rownames(cc)) && unique.names)
    rownames(cc) <- make.unique(
      unlist(lapply(rownames(cc),
                    function(x) toString(x, width=width)))
    )
  if (!std.error) cc <- cc[,-2, drop=FALSE]
  if (!p.value) cc[, -ncol(cc), drop=FALSE]

  sep.pos <- c()
  if (missing(sep.which) && !is.null(x$model.index)) {
    sep.which <- unlist(lapply(x$model.index,
                               function(x)
                                 tail(x,1)))[-length(x$model.index)]
  }
  if (missing(sep.which)) sep.which <- NULL

  if (!is.null(sep.which)) {
    sep0 <- 0%in%sep.which
    if (sep0)
      sep.which <- setdiff(sep.which, 0)
    cc0 <- c()
    sep.which <- c(0, sep.which,nrow(cc))
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
      sep.pos <- c(1,sep.pos+1)
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
    cat("\n", x$compare$method[3],"\n")
    cat(paste(" ", x$compare$method[-(1:3)], collapse="\n"), "\n")
    if (length(x$compare$method)>4) {
      out <- character()
      out <- with(x$compare, c(out, paste(names(statistic),
                                          "=", format(round(statistic, 4)))))
      out <- with(x$compare, c(out, paste(names(parameter),
                                          "=", format(round(parameter,3)))))
      fp  <- with(x$compare, format.pval(p.value, digits = digits))
      out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "<")
                                       fp else paste("=", fp)))
      cat(" ",strwrap(paste(out, collapse = ", ")), sep = "\n")
    }
  }
}

##' @export
vcov.estimate <- function(object, list=FALSE, ...) {
  res <- object$vcov
  nn <- names(coef(object, ...))
  if (list && !is.null(object$model.index)) {
    return(lapply(object$model.index, function(x) object$vcov[x,x]))
  }
  dimnames(res) <- list(nn, nn)
  res
}

##' @export
coef.estimate <- function(object,mat=FALSE,list=FALSE,messages=lava.options()$messages,...) {
  if (mat) return(object$coefmat)
  if (messages>0 && !is.null(object$back.transform)) message("Note: estimates on original scale (before 'back.transform')")
  if (list && !is.null(object$model.index)) {
    return(lapply(object$model.index, function(x) object$coef[x]))
  }
  object$coef
}

##' @export
summary.estimate <- function(object, ...) {
  p <- coef(object, messages=0)
  test <- estimate(coef=p, vcov=vcov(object, messages=0),
                   contrast=as.list(seq_along(p)), ...)
  object$compare <- test$compare
  object <- object[c("coef", "coefmat", "vcov", "call",
                     "ncluster", "model.index", "compare")]
  class(object) <- "summary.estimate"
  object
}

##' @export
coef.summary.estimate <- function(object, ...) {
  object$coefmat
}

##' @export
parameter.estimate <- function(x, ...) {
  return(x$coefmat)
}

##' @export
print.summary.estimate <- function(x, ...) {
  print.estimate(x, type=2L, ...)
}

##' @export
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

##' @export
model.frame.estimate <- function(formula, ...) {
  NULL
}
