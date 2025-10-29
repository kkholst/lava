##' Monte Carlo simulation
##'
##' Applies a function repeatedly for a specified number of replications or over
##' a list/data.frame with plot and summary methods for summarizing the Monte
##' Carlo experiment. Can be parallelized via the future package (use the
##' future::plan function).
##' @export
##' @param x function or 'sim' object
##' @param R Number of replications or data.frame with parameters
##' @param f Optional function (i.e., if x is a matrix)
##' @param colnames Optional column names
##' @param seed (optional) Seed (needed with cl=TRUE)
##' @param args (optional) list of named arguments passed to (mc)mapply
##' @param iter If TRUE the iteration number is passed as first argument to
##'   (mc)mapply
##' @param mc.cores Optional number of cores. Will use parallel::mcmapply
##'   instead of future
##' @param progressr.message Optional message for the progressr progress-bar
##' @param ... Additional arguments to future.apply::future_mapply
##' @aliases sim sim.default as.sim
##' @seealso summary.sim plot.sim print.sim
##' @details To parallelize the calculation use the future::plan function (e.g.,
##'   future::plan(multisession()) to distribute the calculations over the R
##'   replications on all available cores). The output is controlled via the
##'   progressr package (e.g., progressr::handlers(global=TRUE) to enable
##'   progress information).
##' @examples
##' m <- lvm(y~x+e)
##' distribution(m,~y) <- 0
##' distribution(m,~x) <- uniform.lvm(a=-1.1,b=1.1)
##' transform(m,e~x) <- function(x) (1*x^4)*rnorm(length(x),sd=1)
##'
##' onerun <- function(iter=NULL,...,n=2e3,b0=1,idx=2) {
##'     d <- sim(m,n,p=c("y~x"=b0))
##'     l <- lm(y~x,d)
##'     res <- c(coef(summary(l))[idx,1:2],
##'              confint(l)[idx,],
##'              estimate(l,only.coef=TRUE)[idx,2:4])
##'     names(res) <- c("Estimate","Model.se","Model.lo","Model.hi",
##'                     "Sandwich.se","Sandwich.lo","Sandwich.hi")
##'     res
##' }
##' val <- sim(onerun,R=10,b0=1)
##' val
##'
##' val <- sim(val,R=40,b0=1) ## append results
##' summary(val,estimate=c(1,1),confint=c(3,4,6,7),true=c(1,1))
##'
##' summary(val,estimate=c(1,1),se=c(2,5),names=c("Model","Sandwich"))
##' summary(val,estimate=c(1,1),se=c(2,5),true=c(1,1),
##'         names=c("Model","Sandwich"),confint=TRUE)
##'
##' if (interactive()) {
##'     plot(val,estimate=1,c(2,5),true=1,
##'          names=c("Model","Sandwich"),polygon=FALSE)
##'     plot(val,estimate=c(1,1),se=c(2,5),main=NULL,
##'          true=c(1,1),names=c("Model","Sandwich"),
##'          line.lwd=1,col=c("gray20","gray60"),
##'          rug=FALSE)
##'     plot(val,estimate=c(1,1),se=c(2,5),true=c(1,1),
##'          names=c("Model","Sandwich"))
##' }
##'
##' f <- function(a=1, b=1) {
##'   rep(a*b, 5)
##' }
##' R <- Expand(a=1:3, b=1:3)
##' sim(f, R)
##' sim(function(a,b) f(a,b), 3, args=c(a=5,b=5))
##' sim(function(iter=1,a=5,b=5) iter*f(a,b), iter=TRUE, R=5)
sim.default <- function(x = NULL, R = 100, f = NULL, colnames = NULL,
                        seed = NULL, args = list(),
                        iter = FALSE, mc.cores,
                        progressr.message = NULL,
                        ...) {
  stm <- proc.time()
  oldtm <- rep(0, 5)
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    runif(1)
  }
  if (is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  olddata <- NULL
  dots <- list(...)
  mycall <- match.call(expand.dots = FALSE)
  if (inherits(x, c("data.frame", "matrix"))) olddata <- x
  if (inherits(x, "sim")) {
    oldtm <- attr(x, "time")
    oldcall <- attr(x, "call")
    x <- attr(x, "f")
    if (!is.null(f)) x <- f
    ex <- oldcall[["..."]]
    for (nn in setdiff(names(ex), names(dots))) {
      dots[[nn]] <- ex[[nn]]
      val <- list(ex[[nn]])
      names(val) <- nn
      mycall[["..."]] <- c(mycall[["..."]], list(val))
    }
  } else {
    if (!is.null(f)) x <- f
    if (!is.function(x)) stop("Expected a function or 'sim' object.")
  }
  if (is.null(x)) stop("Must give new function argument 'f'.")
  res <- val <- NULL
  on.exit({
    if (is.null(colnames) && !is.null(val)) {
      if (is.matrix(val[[1]])) {
        colnames <- base::colnames(val[[1]])
      } else {
        colnames <- names(val[[1]])
      }
    }
    base::colnames(res) <- colnames
    if (!is.null(olddata)) res <- rbind(olddata, res)
    attr(res, "call") <- mycall
    attr(res, "f") <- x
    cls <- ifelse(is.data.frame(res), "data.frame", "matrix")
    class(res) <- c("sim", cls)
    attr(res, "time") <- proc.time() - stm + oldtm
    return(res)
  })
  if (inherits(R, c("matrix", "data.frame")) || length(R) > 1) {
    if (is.list(R)) {
      parval <- R
      ## list of parameters for each iteration
    } else if (inherits(R, c("matrix", "data.frame"))) {
      parval <- as.data.frame(R)
      names(parval) <- colnames(R)
      R <- NROW(parval)
    }
  } else {
    parval <- as.data.frame(1:R)
    names(parval) <- NULL
  }

  repl <- NROW(parval)
  pb <- progressr::progressor(steps = repl)
  robx <- function(iter__, ...) {
    if (!is.null(progressr.message)) {
      pb(message = progressr.message(...))
    } else {
      pb()
    }
    tryCatch(x(...), error = function(e) NA)
  }
  if (iter || !is.data.frame(parval)) {
    formals(robx)[[1]] <- NULL
  }

  if (is.data.frame(parval)) {
    pp <- c(
      as.list(parval), dots,
      list(FUN = robx, SIMPLIFY = FALSE, MoreArgs = as.list(args))
    )
  } else { ## parameters as a list
    for (i in seq_along(parval)) {
      parval[i] <- c(parval[i], dots, args)
    }
    pp <- c(
      list(parval),
      list(FUN = robx, SIMPLIFY = FALSE)
    )
  }
  if (is.null(pp$future.seed)) {
    pp$future.seed <- TRUE
  }
  if (!missing(mc.cores)) {
    pp$future.seed <- NULL
    pp$mc.cores <- mc.cores
  }

  if (!missing(mc.cores)) {
    if (is.null(mc.cores)) {
      val <- do.call(mapply, pp)
    } else {
      val <- do.call(parallel::mcmapply, pp)
    }
  } else {
    val <- do.call(future.apply::future_mapply, pp)
  }
  res <- do.call(rbind, val)
  if (is.null(res)) {
    res <- matrix(NA, ncol=length(val[[1]]), nrow=repl)
  }
  res
}

##' @export
cbind.sim <- function(x, ...) {
  res <- cbind(as.data.frame(x), ...)
  as.sim(res)
}

##' @export
rbind.sim <- function(x, ...) {
  res <- rbind(as.data.frame(x), ...)
  as.sim(res)
}

##' @export
as.vector.sim <- function(x, mode="any") {
  as.vector(x[,,drop=TRUE], mode=mode)
}

##' @export
as.matrix.sim <- function(x, ...) {
  if (inherits(x, "data.frame")) {
    return(as.matrix(as.data.frame(x)))
  }
  class(x) <- "matrix"
  attr(x, "call") <- NULL
  attr(x, "f") <- NULL
  attr(x, "time") <- NULL
  x
}

##' @export
"[.sim" <- function(x, i, j, drop = FALSE) {
  atr <- attributes(x)
  if (!is.null(dim(x))) {
    class(x) <- class(x)[2]
  } else {
    class(x) <- class(x)[-1]
  }
  x <- NextMethod("[", drop=drop)
  atr.keep <- c("call", "time")
  if (missing(j)) atr.keep <- c(atr.keep, "f")
  attributes(x)[atr.keep] <- atr[atr.keep]
  if (!drop) class(x) <- c("sim", class(x))
  return(x)
}

##' @export
"as.sim" <- function (object, name, ...) {
  if (is.vector(object)) {
    cl <- ifelse(inherits(class(object), "data.frame"), "data.frame", "matrix")
    object <- structure(cbind(object), class=c("sim", cl))
    if (!missing(name)) colnames(object) <- name
    return(object)
  }
  structure(object, class=c("sim", class(object)))
}

Time <- function(sec,print=FALSE,...) {
    h <- sec%/%3600
    m0 <- (sec%%3600)
    m <- m0%/%60
    s <- m0%%60
    res <- c(h=h, m=m, s=s)
    if (print) {
        if (h>0) cat(h, "h ", sep="")
        if (m>0) cat(m, "m ", sep="")
        cat(s, "s", sep="")
        return(invisible(res))
    }
    return(res)
}


##' Generic print method
##'
##' Nicer print method for tabular data. Falls back to standard print method for
##' all other data types.
##' @export
##' @param x object to print
##' @param n number of rows to show from top and bottom of tabular data
##' @param digits precision
##' @param ... additional arguments to print method
Print <- function(x, n=5,
                  digits=max(3, getOption("digits")-3), ...) {
    mat <- !is.null(dim(x))
    if (!mat) {
      if (is.vector(x)) {
          x <- cbind(x)
          colnames(x) <- ""
      } else {
        print(x, ...)
        return(invisible(x))
      }
    } 
    if (is.null(rownames(x))) {
        rownames(x) <- seq(nrow(x))
    }
    sep <- rbind("---"=rep('', ncol(x)))
    if (n<1) {
        print(x, quote=FALSE, digits=digits, ...)
    } else {
      if (NROW(x)<=(2*n)) {
        hd <- base::format(x, digits=digits, ...)
        print(hd, quote=FALSE, ...)
      } else {
        hd <- base::format(utils::head(x, n), digits=digits, ...)
        tl <- base::format(utils::tail(x, n), digits=digits, ...)
        print(rbind(base::as.matrix(hd), sep, base::as.matrix(tl)),
              quote=FALSE, ...)
      }
    }
    invisible(x)
}

##' @export
print.sim <- function(x, ...) {
    s <- summary(x, minimal=TRUE, ...)
    attr(x, "f") <- attr(x, "call") <- NULL
    Print(x, ...)
    cat("\n")
    if (nrow(x)>1) print(s, extra=FALSE, ...)
    return(invisible(x))
}

##' @export
print.summary.sim <- function(x,group=list(c("^mean$","^sd$","^se$","^se/sd$","^coverage"),
                                   c("^min$","^[0-9.]+%$","^max$"),
                                   c("^na$","^missing$"),
                                   c("^true$","^bias$","^rmse$")),
                      lower.case=TRUE,
                      na.print="",
                      digits = max(3, getOption("digits") - 2),
                      quote=FALSE,
                      time=TRUE,
                      extra=TRUE,
                      ...) {
    if (extra) {
        cat(attr(x,"n")," replications",sep="")
        if (time && !is.null(attr(x,"time"))) {
            cat("\t\t\t\t\tTime: ")
            Time(attr(x,"time")["elapsed"],print=TRUE)
        }
        cat("\n\n")
    }

    nn <- rownames(x)
    if (lower.case)  nn <- tolower(nn)
    gg <- lapply(group,
                 function(x) unlist(lapply(x,function(v) grep(v,nn))))
    gg <- c(gg,list(setdiff(seq_along(nn),unlist(gg))))

    x0 <- c()
    ng <- length(gg)
    for (i in seq(ng)) {
        x0 <- rbind(x0, x[gg[[i]],,drop=FALSE],
        { if(i<ng && length(gg[[i+1]])>0) NA})
    }

    print(structure(x0,class="matrix")[,,drop=FALSE],digits=digits,quote=quote,na.print=na.print,...)
    if (extra) cat("\n")
    invisible(x)
}


##' Summary method for 'sim' objects
##'
##' Summary method for 'sim' objects
##' @export
##' @export summary.sim
##' @param object sim object
##' @param estimate (optional) columns with estimates
##' @param se (optional) columns with standard error estimates
##' @param confint (optional) list of pairs of columns with confidence limits
##' @param true (optional) vector of true parameter values
##' @param fun (optional) summary function
##' @param names (optional) names of estimates
##' @param unique.names if TRUE, unique.names will be applied to column names
##' @param minimal if TRUE, minimal summary will be returned
##' @param level confidence level (0.95)
##' @param quantiles quantiles (0,0.025,0.5,0.975,1)
##' @param ... additional levels to lower-level functions
summary.sim <- function(object,estimate=NULL,se=NULL,
                confint=!is.null(se)&&!is.null(true),true=NULL,
                fun,names=NULL,unique.names=TRUE,minimal=FALSE,
                level=0.95,quantiles=c(0,.025,0.5,.975,1),...) {
    if (is.list(estimate)) {
        est <- estimate
        if (is.null(names)) names <- base::names(est)
        estimate <- c()
        nse  <- is.null(se)
        ntrue <- is.null(true)
        elen <- unlist(lapply(est,length))
        est <- lapply(est, function(e) c(e, rep(NA,max(elen)-length(e))))
        for (e in est) {            
            estimate <- c(estimate,e[1])
            if (length(e)>1 && nse) se <- c(se,e[2])
            if (length(e)>2 && ntrue) true <- c(true,e[3])
        }
        cl <- match.call()
        cl[c("estimate","se","true","names")] <- list(estimate,se,true,names)
    }
    if (minimal) {
        fun <- function(x,se,confint,...) {
            res <- c(Mean=mean(x,na.rm=TRUE),
                    SD=sd(x,na.rm=TRUE))
            if (!missing(se) && !is.null(se)) {
                res <- c(res, c(SE=mean(se,na.rm=TRUE)))
                res <- c(res, c("SE/SD"=res[["SE"]]/res[["SD"]]))
            }            
            return(res)
        }
    }
    mfun <- function(x,...) {
        res <- c(mean(x,na.rm=TRUE),
                 sd(x,na.rm=TRUE),
                 if (length(quantiles)>0) quantile(x,quantiles,na.rm=TRUE),
                 mean(is.na(x)))
        if (length(quantiles)>0) {
            nq <- paste0(quantiles*100,"%")
            idx <- which(quantiles==1)
            if (length(idx)>0) nq[idx] <- "Max"
            idx <- which(quantiles==0)
            if (length(idx)>0) nq[idx] <- "Min"
        }
        names(res) <- c("Mean","SD",
                        if (length(quantiles)>0) nq,
                        "Missing")
        res
    }
    tm <- attr(object,"time")
    N <- max(length(estimate),length(se),length(true))    
    if (!is.null(estimate)) estimate <- rep(estimate,length.out=N)
    if (!is.null(se)) se <- rep(se,length.out=N)
    if (!is.null(true)) {
        if (is.null(estimate)) N <- ncol(object)
        true <- rep(true,length.out=N)
    }
    
    if (!is.null(estimate) && is.character(estimate)) {
        estimate <- match(estimate,colnames(object))
    }
    if (!missing(fun)) {
        if (!is.null(estimate)) m.est <- object[,estimate,drop=FALSE]
        else m.est <- object
        m.se <- NULL
        if (!is.null(se)) m.se <- object[,se,drop=FALSE]
        m.ci <- NULL
        if (!is.null(confint)) m.ci <- object[,confint,drop=FALSE]
        res <- lapply(seq(ncol(m.est)),
                      function(i,...) fun(m.est[,i,drop=TRUE],se=m.se[,i,drop=TRUE],confint=m.ci[,1:2+(i-1)*2],...,INDEX=i),...)
        res <- matrix(unlist(res),nrow=length(res[[1]]),byrow=FALSE)
        if (is.null(dim(res))) {
            res <- rbind(res)
        }
        if (is.null(rownames(res))) {
            rownames(res) <- names(fun(0,m.se,m.ci,INDEX=1,...))
            if (is.null(rownames(res))) rownames(res) <- rep("",nrow(res))
        }
        if (is.null(colnames(res))) {
            colnames(res) <- colnames(m.est)
        }
        return(structure(res,
                    n=NROW(object),
                    time=tm,
                    class=c("summary.sim","matrix")))
    }

    if (!is.null(estimate)) {
        est <- apply(object[,estimate,drop=FALSE],2,mfun)
    } else {
        est <- apply(object,2,mfun)
    }

    if (!is.null(true)) {
        if (length(true)!=ncol(est)) {
            ##stop("'true' should be of same length as 'estimate'.")
            true <- rep(true,length.out=ncol(estimate))
        }
        est <- rbind(est,
                     rbind(True=true),rbind(Bias=est["Mean",]-true),
                     rbind(RMSE=((est["Mean",]-true)^2+(est["SD",])^2)^.5)
                     )
    }
    if (!is.null(se)) {
        if (is.character(se)) {
            se <- match(se,colnames(object))
        }
        if (length(se)!=ncol(est)) stop("'se' should be of same length as 'estimate'.")
        est <- rbind(est, SE=apply(object[,se,drop=FALSE],2,
                                   function(x) val <- c(mean(x,na.rm=TRUE))))
        est <- rbind(est,"SE/SD"=est["SE",]/est["SD",])
    }
    if (!is.null(confint) && (length(confint)>1 || confint)) {
        if (is.character(confint)) {
            confint <- match(confint,colnames(object))
        }
        if (length(confint)==1 && confint) {
            if (is.null(se)) stop("Supply confidence limits or SE")
            confint <- c()
            pos <- ncol(object)
            for (i in seq_along(estimate)) {
              z <- 1-(1-level)/2
              if (is.na(se[i])) {
                CI <- matrix(NA, nrow=NROW(object), ncol=2)
              } else {
                CI <- cbind(object[,estimate[i]]-qnorm(z)*object[,se[i]],
                            object[,estimate[i]]+qnorm(z)*object[,se[i]])
              }
              colnames(CI) <- NULL
              object <- cbind(object,CI)
              confint <- c(confint,pos+1:2)
              pos <- pos+2
            }
        }
        if (length(confint)!=2*length(estimate)) stop("'confint' should be of length 2*length(estimate).")
        Coverage <- c()
        for (i in seq_along(estimate)) {
            Coverage <- c(Coverage,
                          mean((object[,confint[2*(i-1)+1]]<true[i]) & (object[,confint[2*i]]>true[i]),na.rm=TRUE))
        }
        est <- rbind(est,Coverage=Coverage)
    }

    if (!is.null(names)) {
         if (length(names)<ncol(est)) {
            uest <- unique(estimate)
            names <- names[match(estimate,uest)]
        }
        colnames(est) <- names

    }
    if (unique.names && !is.null(colnames(est))) {
        colnames(est) <- make.unique(colnames(est))
    }
    est[is.nan(est)] <- NA

    return(structure(est,
                     n=NROW(object),
                     time=tm,
                     class=c("summary.sim","matrix")))
}
