###{{{ gaussian

gaussian_method.lvm <- "nlminb2"
`gaussian_objective.lvm` <-
    function(x,p,data,S,mu,n,...) {
        mp <- modelVar(x,p=p,data=data,...)
        C <- mp$C ## Model specific covariance matrix
        xi <- mp$xi ## Model specific mean-vector
        if (!lava.options()$allow.negative.variance && any(diag(mp$P)<0)) return(NaN)

        iC <- Inverse(C,det=TRUE, symmetric = TRUE)
        detC <- attributes(iC)$det
        if (n<2) {
            z <- as.numeric(data-xi)
            val <- log(detC) + tcrossprod(z,crossprod(z,iC))[1]
            return(0.5*val)
        }
        if (!is.null(mu)){
            W <- suppressMessages(crossprod(rbind(mu-xi)))
            T <- S+W
        } else {
            T <- S
        }
        res <- n/2*log(detC) + n/2*tr(T%*%iC) ## Objective function (Full Information ML)
        ## if (any(attr(iC,"lambda")<1e-16)) res <- res-1e2
        return(res)
    }


`gaussian_hessian.lvm` <- function(x,p,n,...) {
    dots <- list(...); dots$weights <- NULL
    do.call("information", c(list(x=x,p=p,n=n),dots))
}

gaussian_gradient.lvm <-  function(x,p,data,S,mu,n,...) {
    dots <- list(...); dots$weights <- NULL
    if (n>2) data <- NULL
    val <- -gaussian_score.lvm(x,p=p,S=S,mu=mu,n=n,data=data,reindex=FALSE,...)
    if (!is.null(nrow(val))) {
        val <- colSums(val)
    }
    val
}

gaussian_score.lvm <- function(x, data, p, S, n, mu=NULL, weights=NULL, debug=FALSE, reindex=FALSE, mean=TRUE, constrain=TRUE, indiv=FALSE,...) {

    if (!is.null(data)) {
        if ((nrow(data)<2 | !is.null(weights))| indiv)
        {
            mp <- modelVar(x,p,data=data[1,])
            iC <- Inverse(mp$C,det=FALSE, symmetric = TRUE)
            MeanPar <- attributes(mp)$meanpar
            D <- with(attributes(mp), deriv.lvm(x, meanpar=MeanPar, p=pars, mom=mp, mu=NULL)) ##, all=length(constrain(x))>0))
            myvars <- (index(x)$manifest)
            data <- data[,myvars,drop=FALSE]
            score <- matrix(ncol=length(p),nrow=NROW(data))
            score0 <- -1/2*as.vector(iC)%*%D$dS
            if (!is.null(weights)) {
                W0 <- diag(nrow=length(myvars))
                widx <- match(colnames(weights),myvars)
            }

            for (i in seq_len(NROW(data))) {
                z <- as.numeric(data[i,])
                u <- z-as.numeric(mp$xi)
                if (!is.null(weights)) {
                    W <- W0; diag(W)[widx] <- as.numeric(weights[i,])
                    score[i,] <-
                        as.numeric(crossprod(u,iC%*%W)%*%D$dxi +
                                   -1/2*(as.vector((iC
                                       - iC %*% crossprod(rbind(u))
                                       %*% iC)%*%W)) %*% D$dS
                                   )
                } else {
                  score[i, ] <-
                    as.numeric(score0 + crossprod(u, iC)%*%D$dxi +
                      1/2*as.vector(iC%*%crossprod(rbind(u))%*%iC)%*%D$dS)
                }
            }; colnames(score) <- names(p)
            return(score)
        }
    }

    ### Here the emperical mean and variance of the population are sufficient statistics:
    if (missing(S)) {
        data0 <- na.omit(data[,manifest(x),drop=FALSE])
        n <- NROW(data0)
        S <- cov(data0)*(n-1)/n
        mu <- colMeans(data0)
    }
    mp <- modelVar(x,p)
    C <- mp$C
    xi <- mp$xi
    iC <- Inverse(C,det=FALSE, symmetric = TRUE)
    Debug("Sufficient stats.",debug)
    if (!is.null(mu) & !is.null(xi)) {
        W <- crossprod(rbind(mu-xi))
        T <- S+W
    } else {
        T <- S
    }
    D <- deriv.lvm(x, meanpar=attributes(mp)$meanpar, mom=mp, p=p, mu=mu, mean=mean)
    vec.iC <- as.vector(iC)
    if (lava.options()$devel) {
        Grad <- numeric(length(p))
        imean <- with(index(x)$parBelongsTo,mean)
        Grad[-imean] <- n/2*crossprod(D$dS[,-imean], as.vector(iC%*%T%*%iC)-vec.iC)
    } else {
        Grad <- n/2*crossprod(D$dS, as.vector(iC%*%T%*%iC)-vec.iC)
    }
    if (!is.null(mu) & !is.null(xi)) {
        if (!(lava.options()$devel)) {
            Grad <- Grad - (n/2*crossprod(D$dT,vec.iC))
        } else {
            Grad[with(index(x)$parBelongsTo,c(mean,reg))] <- Grad[with(index(x)$parBelongsTo,c(mean,reg))] - (n/2*crossprod(D$dT,vec.iC))
        }
    }
    res <- as.numeric(Grad)
    return(rbind(res))
}

###}}} gaussian

###{{{ gaussian variants

## Maximum Likelihood with numerical gradient + hessian
gaussian0_objective.lvm <- gaussian_objective.lvm

gaussian1_objective.lvm <- gaussian_objective.lvm
gaussian1_gradient.lvm <- function(...) gaussian_gradient.lvm(...)
gaussian1_hessian.lvm <- function(x,p,...) {
    myg2 <- function(p1) gaussian_gradient.lvm(x,p=p1,...)
    myg3 <- function(p1) numDeriv::jacobian(myg2,p1)
    return(myg3(p))
    ## myg <- function(p1) gaussian_objective.lvm(x,p=p1,...)
    ## numDeriv::hessian(myg,p)
}

## BHHH
gaussian2_method.lvm <- "NR"
gaussian2_objective.lvm <- gaussian_objective.lvm
gaussian2_gradient.lvm <- gaussian_gradient.lvm
gaussian2_hessian.lvm <- function(x,p,n,data,...) {
    S <- -score(x,p=p,n=n,data=data,indiv=TRUE,...)
    I <- t(S)%*%S
    attributes(I)$grad <- colSums(S)
    return(I)
}

###}}}

###{{{ Weighted

weighted_method.lvm <- "NR"
weighted_gradient.lvm <- function(x,p,data,weights,indiv=FALSE,...) {
    myvars <- index(x)$manifest
    if (NCOL(data)!=length(myvars))
        data <- subset(data,select=myvars)
    score <- matrix(ncol=length(p),nrow=NROW(data))
    myy <- index(x)$endogenous
    myx <- index(x)$exogenous
    mynx <- setdiff(myvars,myx)
    W0 <- diag(nrow=length(myy))
    widx <- match(colnames(weights),myy)
    pp <- modelPar(x,p)
    mp <- moments(x,p=p,conditional=TRUE,data=data[1,])
    iC <- Inverse(mp$C,det=FALSE, symmetric = TRUE)
    v <- matrix(0,ncol=length(vars(x)),nrow=NROW(data))
    colnames(v) <- vars(x)
    for (i in mynx) v[,i] <- mp$v[i]
    for (i in myx) v[,i] <- data[,i]
    xi <- t(mp$G%*%t(v))
    u <- as.matrix(data)[,myy]-xi
    D <- deriv.lvm(x, meanpar=pp$meanpar,
                   p=pp$p, mom=mp, mu=NULL)
    if (NROW(data)==1) {
        W <- W0; diag(W)[widx] <- as.numeric(weights[i,])
        score[i,] <-
            as.numeric(crossprod(u,iC%*%W)%*%D$dxi +
                       -1/2*(as.vector((iC
                           - iC %*% crossprod(rbind(u))
                           %*% iC)%*%W)) %*% D$dS)
        return(-score)
    }
    Gdv <- mp$G%*%D$dv
    for (i in seq_len(NROW(data))) {
        W <- W0; diag(W)[widx] <- as.numeric(weights[i,])
        dxi <-
            (t(as.numeric(v[i,]))%x%diag(nrow=length(myy)))%*%D$dG + Gdv
        score[i,] <- -0.5*as.vector(iC%*%W)%*%D$dS +
            as.numeric(crossprod(u[i,],iC%*%W)%*%dxi +
                       1/2*as.vector(iC%*%crossprod(rbind(u[i,]))%*%iC%*%W)%*%D$dS)
        ## score[i,] <- -0.5*as.vector(iC)%*%D$dS +
        ##   as.numeric(crossprod(u[i,],iC)%*%dxi +
        ##              1/2*as.vector(iC%*%tcrossprod(u[i,])%*%iC)%*%D$dS)

    }
    if (indiv) return(-score)
    colSums(-score)
}
weighted_hessian.lvm <- function(...) {
    S <- weighted_gradient.lvm(...,indiv=TRUE)
    res <- crossprod(S)
    attributes(res)$grad <- colSums(-S)
    res
}


weighted0_method.lvm <- "estimatingfunction"
weighted0_gradient.lvm <- function(...) {
    val <- -gaussian_score.lvm(...)
    colSums(val)
}
weighted0_hessian.lvm <- NULL

weighted2_method.lvm <- "estimatingfunction"
weighted2_gradient.lvm <- function(x,p,data,weights,indiv=FALSE,...) {
    myvars <- index(x)$manifest
    if (NCOL(data)!=length(myvars))
        data <- subset(data,select=myvars)
    score <- matrix(ncol=length(p),nrow=NROW(data))
    myy <- index(x)$endogenous
    myx <- index(x)$exogenous
    mynx <- setdiff(myvars,myx)
    W0 <- diag(nrow=length(myy))
    widx <- match(colnames(weights),myy)
    pp <- modelPar(x,p)
    for (i in seq_len(NROW(data))) {
        z <- as.matrix(data[i,myy])
        mp <- moments(x,p=p,conditional=TRUE,data=data[i,])
        u <- as.numeric(z-mp$xi[,1])
        iC <- Inverse(mp$C,det=FALSE, symmetric = TRUE)
        D <- deriv.lvm(x, meanpar=pp$meanpar,
                       p=pp$p, mom=mp, mu=NULL)
        W <- W0; diag(W)[widx] <- as.numeric(weights[i,])
        score[i,] <- -0.5*as.vector(iC%*%W)%*%D$dS +
            as.numeric(crossprod(u,iC%*%W)%*%D$dxi +
                       1/2*as.vector(iC%*%crossprod(rbind(u))%*%iC%*%W)%*%D$dS)
    }
    if (indiv) return(-score)
    colSums(-score)
}
weighted2_hessian.lvm <- NULL

###}}} Weighted

###{{{ Simple

`Simple_hessian.lvm` <- function(p,...) {
    matrix(NA, ncol=length(p), nrow=length(p))
}
Simple_gradient.lvm <- function(x,p,...) {
    naiveGrad(function(pp) Simple_objective.lvm(x,pp,...), p)
}
`Simple_objective.lvm` <-
    function(x, p=p, S=S, n=n, ...) {
        m. <- moments(x,p)
        C <- m.$C
        npar <- m.$npar
        detC <- det(C)
        iC <- Inverse(C, symmetric = TRUE)
        if (detC<0 | inherits(iC, "try-error"))
            return(.Machine$double.xmax)
        res <- n/2*(log(detC) + tr(S%*%iC) - log(det(S)) - npar)
        res
    }

###}}} ObjectiveSimple
