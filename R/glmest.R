
glm_estimate_hook <- function(x,estimator,...) {
    yy <- c()
    if (length(estimator)>0 && estimator=="glm") {
        for (y in endogenous(x)) {
            fam <- attributes(distribution(x)[[y]])$family
            if (is.null(fam)) fam <- stats::gaussian()
            if (!(tolower(fam$family)%in%
                  c("gaussian","gamma","inverse.gaussian","weibull"))) {
                yy <- c(yy,y)
            }
        }
        if (length(yy)>0) covariance(x,yy) <- 1
    }
    return(c(list(x=x,estimator=estimator,...)))
}

GLMest <- function(m,data,control=list(),...) {
    yvar <- endogenous(m)
    res <- c()
    count <- 0
    V <- NULL
    mymsg <- c()
    ics <- c()
    breads <- c()

    et <- eventTime(m)
    yvar.et <- rep(NA,length(yvar))
    names(yvar.et) <- yvar
    if (!is.null(et)) {
        for (i in seq_along(et)) {
            ## if (!survival::is.Surv(data[,et[[i]]$names[1]]))
            ##         data[,et[[i]]$names[1]] <- with(et[[i]],
            ##         survival::Surv(data[,names[1]],data[,names[2]]))
            yvar <- setdiff(yvar,c(et[[i]]$latentTimes[-1],et[[i]]$names))
            yvar.et[et[[i]]$latentTimes[1]] <- et[[i]]$names[1]
        }
    }

    for (y in yvar) {
        count <- count+1
        xx <- parents(m, y)
        fam <- attributes(distribution(m)[[y]])$family
        if (is.null(fam)) fam <- stats::gaussian()
        if (!is.null(fam$link)) {
            mymsg <- c(mymsg, with(fam, paste0(family,"(",link,")")))
        } else {
            mymsg <- c(mymsg, with(fam, paste0(family)))
        }
        if (length(xx)==0) xx <- 1
        nn0 <- paste(y,xx,sep=lava.options()$symbol[1])
        y0 <- y
        ## isEventTime <- !is.na(yvar.et[y])
        ## if (isEventTime) {
        ##     y <- yvar.et[y]
        ## }
        # nn0 <- paste(y,xx,sep=lava.options()$symbol[1])

        f <- as.formula(paste0(y, "~",
          paste(xx, collapse = "+")
        ))
        isSurv <- inherits(data[1, y], "Surv")
        if (isSurv) {
          g <- survival::survreg(f,data=data,dist=fam$family,...)
        } else {
          g <- glm(f,family=fam,data=data,...)
        }
        p <- pars(g)
        ii <- IC(g)
        V0 <- attr(ii, "bread")
        ics <- cbind(ics,ii)
        y <- y0
        names(p)[1] <- y
        if (length(p)>1) {
          xx0 <- setdiff(xx, 1)
          if (length(xx0) > 0) {
            nn0 <- paste(y, xx0, sep = lava.options()$symbol[1])
            names(p)[seq_along(nn0) + 1] <- nn0
          }
          if (length(p)>length(xx0)+1) names(p)[length(p)] <- paste(y,y,sep=lava.options()$symbol[2])
        }
        if (tolower(fam$family)%in%c("gaussian","gamma","inverse.gaussian") && !isSurv) {
            ics <- cbind(ics,0)
            null <- matrix(0); dimnames(null) <- list("scale","scale")
            V0 <- blockdiag(V0,null,pad=0)
        }
        breads <- c(breads, list(V0))
        res <- c(res, list(p));
    }
    coefs <- unlist(res)
    idx <- na.omit(match(coef(m),names(coefs)))
    coefs <- coefs[idx]
    V <- var_ic(ics[, idx])
    mymsg <- noquote(cbind(mymsg))
    colnames(mymsg) <- "Family(Link)"; rownames(mymsg) <- paste(yvar,":")
    list(estimate=coefs, vcov=V, breads=breads, IC=ics[,idx],
         summary.message=function(...)  { mymsg }, dispname="Dispersion:")
}

GLMscore <- function(x,p,data,indiv=TRUE,logLik=FALSE,...) {
    yvar <- endogenous(x)
    S <- pnames <- c()
    count <- 0
    breads <- c()
    L <- 0
    for (y in yvar) {
        count <- count+1
        xx <- parents(x,y)
        pname <- c(y,paste0(y,sep=lava.options()$symbol[1],xx),paste(y,y,sep=lava.options()$symbol[2]))
        pidx <- na.omit(match(pname,coef(x)))
        fam <- attributes(distribution(x)[[y]])$family
        if (is.null(fam)) fam <- stats::gaussian()
        if (length(xx)==0) xx <- 1
        f <- as.formula(paste0(y,"~",paste(xx,collapse="+")))
        isSurv <- inherits(data[1,y],"Surv")
        if (inherits(data[,y],"Surv")) {
            g <- survival::survreg(f,data=data,dist=fam$family)
        } else {
            g <- glm(f,family=fam,data=data)
        }
        p0 <- p[pidx]
        if (!isSurv) L0 <- logL.glm(g,p=p0,indiv=TRUE,...)
        if (tolower(fam$family)%in%c("gaussian","gamma","inverse.gaussian") && !isSurv) {
            p0 <- p0[-length(p0)]
            S0 <- score(g,p=p0,indiv=TRUE,pearson=TRUE,...)
            V0 <- attr(S0,"bread")
            r <- attr(S0,"pearson")
            ## dispersion <- mean(r^2)
            S0 <- cbind(S0,scale=0)
            null <- matrix(0); dimnames(null) <- list("scale","scale")
            V0 <- blockdiag(V0,null,pad=0)
        } else {
            S0 <- score(g,p=p0,indiv=TRUE,...)
            if (isSurv) L0 <- attr(S0,"logLik")
            V0 <- attr(S0,"bread")
        }
        L <- L+sum(L0)
        breads <- c(breads,list(V0))
        S <- c(S,list(S0))
        pnames <- c(pnames, list(pname));
    }
    coefs <- unlist(pnames)
    idx <- na.omit(match(coefs,coef(x)))
    idx <- order(idx)
    V <- Reduce(blockdiag,breads)[idx,idx]
    S1 <- Reduce(cbind,S)[,idx,drop=FALSE]
    colnames(S1) <- coef(x)
    attributes(S1)$bread <- V
    attributes(S1)$logLik <- structure(L,nobs=nrow(data),nall=nrow(data),df=length(p),class="logLik")
    if (!indiv) S1 <- colSums(S1)
    return(S1)
}


##' @export
score.lm <- function(x, p=coef(x), data, indiv=FALSE,
              y, X, offset=NULL, weights=NULL, dispersion=TRUE, ...) {
    if (missing(data)) {
        X <- model.matrix(x)
        y <- model.frame(x)[,1]
    } else {
        X <- model.matrix(formula(x),data=data)
        y <- model.frame(formula(x),data=data)[,1]
    }
    if(any(is.na(p))) warning("Over-parameterized model")
    Xbeta <- X%*%p
    if (is.null(offset)) offset <- x$offset
    if (!is.null(offset)) Xbeta <- Xbeta+offset
    r <- y-Xbeta
    if (is.null(weights)) weights <- x$weights
    if (!is.null(weights)) {
        sigma2 <- sum(r^2*weights)/(length(r)-length(p))
        r <- r*weights
    } else {
        sigma2 <- sum(r^2)/(length(r)-length(p))
    }
    if (!dispersion) sigma2 <- 1
    ##sigma2 <- suppressWarnings(summary(x)$sigma^2)
    A <- as.vector(r)/sigma2
    S <- apply(X,2,function(x) x*A)
    if (!indiv) return(colSums(S))
    suppressWarnings(attributes(S)$bread <- vcov(x)*NROW(S))
    return(S)
}

##' @export
score.glm <- function(x,p=coef(x),data,indiv=FALSE,pearson=FALSE,
               y,X,link,dispersion,offset=NULL,weights=NULL,...) {

    if (inherits(x,"glm")) {
        link <- family(x)
        if (missing(data)) {
            X <- model.matrix(x)
            y <- model.frame(x)[,1]
        } else {
            X <- model.matrix(formula(x),data=data)
            y <- model.frame(formula(x),data=data)[,1]
        }
        offset <- x$offset
    } else {
        if (missing(link)) stop("Family needed")
        if (missing(data)) stop("data needed")
        X <- model.matrix(formula(x),data=data)
        y <- model.frame(formula(x),data=data)[,1]
    }
    if (is.character(y) || is.factor(y)) {
        y <- as.numeric(as.factor(y)) - 1
    }
    ## g <- link$linkfun
    ginv <- link$linkinv
    dginv <- link$mu.eta ## D[linkinv]
    ##dg <- function(x) 1/dginv(g(x)) ## Dh^-1 = 1/(h'(h^-1(x)))
    if (inherits(x, "negbin")) {
        Dcanlink <- function(x) 1/x
    } else {
        canonf <- do.call(link$family,list())
        ## caninvlink <- canonf$linkinv
        canlink <- canonf$linkfun
        Dcaninvlink <- canonf$mu.eta
        Dcanlink <- function(x) 1 / Dcaninvlink(canlink(x))
    }
    ##gmu <- function(x) g(caninvlink(x))
    ##invgmu <- function(z) canlink(ginv(z))
    h <- function(z) Dcanlink(ginv(z))*dginv(z)
    if(any(is.na(p))) stop("Over-parameterized model")
    Xbeta <- X%*%p
    if (!is.null(offset)) Xbeta <- Xbeta+offset
    if (missing(data) && !is.null(x$offset) && is.null(offset) ) Xbeta <- Xbeta+x$offset
    pi <- ginv(Xbeta)
    r <- y-pi
    if (!is.null(x$prior.weights) || !is.null(weights)) {
        if (is.null(weights)) weights <- x$prior.weights
    } else {
        weights <- !is.na(r)
    }
    a.phi <- 1
    r <- r*weights
    rpearson <- as.vector(r)/link$variance(pi)^.5
    if (length(p) > length(coef(x))) {
        a.phi <- p[length(coef(x)) + 1]
    } else { ## if (tolower(family(x)$family)%in%c("gaussian","gamma","inverse.gaussian"))
        suppressWarnings(disp <- summary(x)$dispersion)
        if (!is.null(disp)) a.phi <- disp
        if (inherits(x, "negbin")) {
            a.phi <- link$variance(pi) / pi
        }
        ## a.phi <- sum(rpearson^2)*x$df.residual/x$df.residual^2
    }
    A <- as.vector(h(Xbeta) * r) / a.phi
    S <- apply(X, 2, function(x) x * A)
    if (!indiv) return(colSums(S))
    if (pearson) attr(S,"pearson") <- rpearson
    suppressWarnings(attributes(S)$bread <- vcov(x)*NROW(S))
    if (x$family$family=="quasi" && x$family$link=="identity" && x$family$varfun=="constant")
        attributes(S)$bread <- -Inverse(information.glm(x)*NROW(S))
    return(S)
}

##' @export
pars.glm <- function(x,...) {
    if (tolower(family(x)$family)%in%c("gaussian","gamma","inverse.gaussian")) {
        res <- c(coef(x),suppressWarnings(summary(x)$dispersion))
        names(res)[length(res)] <- "Dispersion"
        return(res)
    }
    return(coef(x))
}

logL.glm <- function(x,p=pars.glm(x),data,indiv=FALSE,...) {
    if (!missing(data)) {
        x <- update(x,data=data,...)
    }
    f <- family(x)
    ginv <- f$linkinv
    X <- model.matrix(x)
    n <- nrow(X)
    ##disp <- 1;
    p0 <- p
    if (tolower(family(x)$family)%in%c("gaussian","gamma","inverse.gaussian")) {
        if (length(p)==ncol(X)) {
            ##disp <- suppressWarnings((summary(x)$dispersion))
        } else {
            ##disp <- tail(p,1)
            p0 <- p[-length(p)]
        }
    }
    if(any(is.na(p))) {
        warning("Over-parametrized model")
    }
    Xbeta <- X%*%p0
    if (!is.null(x$offset)) Xbeta <- Xbeta+x$offset
    y <- model.frame(x)[,1]
    mu <- ginv(Xbeta)
    w <- x$prior.weights
    dev <-  f$dev.resids(y,mu,w)
    if (indiv) {

    }
    loglik <- length(p)-(f$aic(y,n,mu,w,sum(dev))/2+x$rank)
    structure(loglik,nobs=n,df=length(p),class="logLik")
}

##' @export
IC.glm <- function(x,...) {
    IC.default(x,...)
}

hessian.glm <- function(x,p=coef(x),...) {
    numDeriv::jacobian(function(theta) score.glm(x,p=theta,indiv=FALSE,...),p)
}

##' @export
information.glm <- function(x,...) hessian.glm(x,...)

robustvar <- function(x,id=NULL,...) {
    U <- score(x,indiv=TRUE)
    II <- unique(id)
    K <- length(II)
    J <- 0
    if (is.null(id)) {
        J <- crossprod(U)
    } else {
        for (ii in II) {
            J <- J+tcrossprod(colSums(U[which(id==ii),,drop=FALSE]))
        }
        J <- K/(K-1)*J
    }
    iI <- vcov(x)
    V <- iI%*%J%*%iI
    return(V)
}

glm_logLik.lvm <- function(object,...) {
    attr(GLMscore(object,...),"logLik")
}

glm_method.lvm <- NULL
glm_objective.lvm <- function(x,p,data,...) {
    GLMest(x,data,...)
}
glm_gradient.lvm <- function(x,p,data,...) {
    -GLMscore(x,p,data,...)
}

glm_variance.lvm <- function(x,p,data,opt,...) {
    opt$vcov
}
