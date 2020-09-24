###{{{ mixture

#' Estimate mixture latent variable model.
#'
#' Estimate mixture latent variable model
#'
#' Estimate parameters in a mixture of latent variable models via the EM
#' algorithm.
#'
#' The performance of the EM algorithm can be tuned via the \code{control}
#' argument, a list where a subset of the following members can be altered:
#'
#' \describe{ \item{start}{Optional starting values} \item{nstart}{Evaluate
#' \code{nstart} different starting values and run the EM-algorithm on the
#' parameters with largest likelihood} \item{tol}{Convergence tolerance of the
#' EM-algorithm.  The algorithm is stopped when the absolute change in
#' likelihood and parameter (2-norm) between successive iterations is less than
#' \code{tol}} \item{iter.max}{Maximum number of iterations of the
#' EM-algorithm} \item{gamma}{Scale-down (i.e. number between 0 and 1) of the
#' step-size of the Newton-Raphson algorithm in the M-step} \item{trace}{Trace
#' information on the EM-algorithm is printed on every \code{trace}th
#' iteration} }
#'
#' Note that the algorithm can be aborted any time (C-c) and still be saved
#' (via on.exit call).
#'
 #' @param x List of \code{lvm} objects. If only a single \code{lvm} object is
#' given, then a \code{k}-mixture of this model is fitted (free parameters
#' varying between mixture components).
#' @param data \code{data.frame}
#' @param k Number of mixture components
#' @param control Optimization parameters (see details)
#' #type Type of EM algorithm (standard, classification, stochastic)
#' @param vcov of asymptotic covariance matrix (NULL to omit)
#' @param names If TRUE returns the names of the parameters (for defining starting values)
#' @param ... Additional arguments parsed to lower-level functions
#' @author Klaus K. Holst
#' @seealso \code{mvnmix}
#' @keywords models regression
#' @examples
#'
#' \donttest{
#' m0 <- lvm(list(y~x+z,x~z))
#' distribution(m0,~z) <- binomial.lvm()
#' d <- sim(m0,2000,p=c("y~z"=2,"y~x"=1),seed=1)
#'
#' ## unmeasured confounder example
#' m <- baptize(lvm(y~x, x~1));
#' intercept(m,~x+y) <- NA
#'
#' if (requireNamespace('mets', quietly=TRUE)) {
#'   set.seed(42)
#'   M <- mixture(m,k=2,data=d,control=list(trace=1,tol=1e-6))
#'   summary(M)
#'   lm(y~x,d)
#'   estimate(M,"y~x")
#'   ## True slope := 1
#' }
#' }
#'
#' @export mixture
mixture <- function(x, data, k=length(x),
             control=list(),
             vcov="observed",
             names=FALSE,
             ...) {
    MODEL <- "normal"
    ##type=c("standard","CEM","SEM"),
    ## type <- tolower(type[1])
    ## if (type[1]!="standard") {
    ##     return(mixture0(x,data=data,k=k,control=control,type=type,...))
    ## }

    optim <- list(
        rerun=TRUE,
        ## EM (accelerated):
        K=1, ## K=1, first order approx, slower but some times more stable
        square=TRUE,
        step.max0=1,
        step.min0=1,
        mstep=4,
        kr=1,
        objfn.inc=2,
        keep.objfval=TRUE,
        convtype= "parameter",
        ## convtype = "objfn",
        maxiter=1500,
        tol=1e-5,
        trace=0,
        ## Starting values:
        start=NULL,
        startbounds=c(-2,2),
        startmean=FALSE,
        nstart=1,
        prob=NULL,
        ## Newton raphson:
        delta=1e-2,
        constrain=TRUE,
        stopc=2,
        lbound=1e-9,
        stabil=TRUE,
        gamma=0.5,
        gamma2=1,
        newton=10,
        lambda=0 # Stabilizing factor (avoid singularities of I)
    )

    if (!missing(control))
        optim[base::names(control)] <- control
    if ("iter.max"%in%base::names(optim)) optim$maxiter <- optim$iter.max

    if (k==1) {
        if (is.list(x))
            res <- estimate(x[[1]],data,...)
        else
            res <- estimate(x,data,...)
        return(res)
    }

    start0 <- NULL

    xx <- x
    if (inherits(x,"lvm")) {
        xx <- rep(list(x), k)
    }
    mg <- multigroup(xx, rep(list(data),k), fix=FALSE)
    ppos <- parpos(mg)
    parname <- attr(ppos,"name")
    naparname <- which(is.na(parname))
    parname[naparname]  <- mg$name[naparname]
    if (names) {
        return(parname)
    }

    if (class(x)[1]=="lvm") {
        index(x) <- reindex(x,zeroones=TRUE,deriv=TRUE)
        if ((is.null(optim$start) || length(optim$start)<length(parname))) {
            start0 <- tryCatch(estimate(x,data,quick=TRUE),
                               error=function(...) NULL)
        }
    }

    ## Bounds on variance parameters
    Npar <- with(mg, npar+npar.mean)
    ParPos <- modelPar(mg,1:Npar)$p
    lower <- rep(-Inf, mg$npar);
    offdiagpos <- c()
    varpos <- c()
    for (i in 1:k) {
        vpos <- sapply(mg$parlist[[i]][variances(mg$lvm[[i]])], function(y) as.numeric(substr(y,2,nchar(y))))
        offpos <- sapply(mg$parlist[[i]][offdiags(mg$lvm[[i]])], function(y) as.numeric(substr(y,2,nchar(y))))
        varpos <- c(varpos, vpos)
        offdiagpos <- c(offdiagpos,offpos)
        if (length(vpos)>0)
            lower[vpos] <- optim$lbound  ## Setup optimization constraints
    }
    lower <- c(rep(-Inf,mg$npar.mean), lower)
    constrained <- which(is.finite(lower))
    if (!any(constrained)) optim$constrain <- FALSE

    mymodel <- list(multigroup=mg,k=k,data=data,parpos=ParPos); class(mymodel) <- "lvm.mixture"

    selfstarter <- FALSE
    if (any(is.na(optim$start)) || length(optim$start)<length(ppos)) {
        constrLogLikS <- function(p) {
            if (optim$constrain) {
                p[constrained] <- exp(p[constrained])
            }
            -logLik(mymodel,c(p=p,rep(1/k,k-1)))
        }

        start <- rep(NA,Npar)
        if (!is.null(start0)) {
            start[ParPos[[1]]] <- start0[seq_along(ParPos[[1]])]
            ii <- which(is.na(start))
            start[ii] <- runif(length(ii),optim$startbounds[1],optim$startbounds[2])
        } else {
            selfstarter <- TRUE
            start <- runif(Npar,optim$startbounds[1],optim$startbounds[2]);
            if (length(offdiagpos)>0)
                start[mg$npar.mean + offdiagpos] <- 0
            if (optim$nstart>1) {
                myll <- constrLogLikS(start)
                for (i in 1:optim$nstart) {
                    newstart <- runif(Npar,optim$startbounds[1],optim$startbounds[2]);
                    newmyll <- constrLogLikS(newstart)
                    if (newmyll<myll) {
                        start <- newstart
                    }
                }
            }
        }
        if (!is.null(optim$start)) {
            iistart <- match(names(optim$start), parname)
            start[iistart] <- optim$start
        }
        optim$start <- start
    }
    if (length(optim$start)>Npar) {
        optim$prob <- optim$start[Npar+seq_len(k-1)]
        optim$start <- optim$start[seq_len(Npar)]
    }

    if (is.null(optim$prob))
        optim$prob <- rep(1/k,k-1)
    thetacur <- optim$start

    probcur <- with(optim, c(prob,1-sum(prob)))

    if (optim$constrain & !selfstarter) {
        thetacur[constrained] <- log(thetacur[constrained])
    }

    PosteriorProb <- function(pp,priorprob,constrain=FALSE) {
        if (!is.list(pp)) {
            if (constrain) {
                pp[constrained] <- exp(pp[constrained])
            }
            if (missing(priorprob)) priorprob <- pp[seq(Npar+1,length(pp))]
            pp <- lapply(ParPos,function(x) pp[x])
        }
        priorprob <- pmax(priorprob,1e-16)
        priorprob <- priorprob/sum(priorprob)
        k <- length(pp)
        logff <- sapply(seq(k), function(j) (logLik(mg$lvm[[j]],p=pp[[j]],data=data,indiv=TRUE,model=MODEL)))
        logplogff <- t(apply(logff,1, function(z) z+log(priorprob)))
        ## Log-sum-exp (see e.g. NR)
        zmax <- apply(logplogff,1,max)
        logsumpff <- log(rowSums(exp(logplogff-zmax)))+zmax
        gamma <- exp(apply(logplogff,2,function(y) y - logsumpff)) ## Posterior class probabilities
        return(gamma)
    }

    negLogLik <- function(p) {
        if (optim$constrain) {
            p[constrained] <- exp(p[constrained])
        }
        myp <- lapply(ParPos,function(x) p[x])
        K <- length(myp)
        prob <- p[seq(Npar+1,Npar+K-1)]; prob <- c(prob,1-sum(prob))
        logff <- sapply(1:length(myp), function(j) (logLik(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE,model=MODEL)))
        ## logff <- sapply(1:length(myp),
        ##              function(j) -normal_objective.lvm(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE))
        logplogff <- t(apply(logff,1, function(y) y+log(prob)))
        zmax <- apply(logplogff,1,max)
        logsumpff <- log(rowSums(exp(logplogff-zmax)))+zmax
        loglik <- sum(logsumpff)
        return(-loglik)
    }


    ObjEstep <- function(p,gamma,pr) {
        if (optim$constrain) {
            p[constrained] <- exp(p[constrained])
        }
        myp <- lapply(ParPos,function(x) p[x])
        loglik = 0;
        for (j in seq_along(myp))
            loglik = loglik + gamma[,j]*(logLik(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE,model=MODEL))
        ## logff <- sapply(1:length(myp), function(j) gamma[,j]*(logLik(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE,model=MODEL)))
        ## zmax <- apply(logff,1,max)
        ## ffz <- apply(logff,2,function(x) exp(x-zmax))
        ## logsff <- log(rowSums(ffz))+zmax
        ## loglik <- sum(logsff)
        return(-sum(loglik))
    }

    GradEstep <- function(p,gamma,pr) {
        if (optim$constrain) {
            p[constrained] <- exp(p[constrained])
        }
        myp <- lapply(ParPos,function(x) p[x])
        ## logff <- sapply(1:length(myp), function(j) log(gamma[,j])+(logLik(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE,model=MODEL)))
        ## Exp-sum normalization: exp(xi)/sum(exp(xi)) = exp(xi-b)*exp(b)/sum(exp(xi-b)exp(b))
        ## = exp(xi-b)/sum(exp(xi-b)),  b=max(xi)
        ## zmax <- apply(logff,1,max)
        ## ffz <- apply(logff,2,function(x) exp(x-zmax))
        ## sffz <- rowSums(ffz)
        D <- lapply(1:length(myp), function(j) {
            ## K <- ffz[,j]/sffz
            val <- score(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE,model=MODEL)
            apply(val,2,function(x) x*gamma[,j])
        })
        D0 <- matrix(0,nrow(data),length(p))
        for (j in 1:k) D0[,ParPos[[j]]] <- D0[,ParPos[[j]]]+D[[j]]
        S <- colSums(D0)
        if (optim$constrain) {
            S[constrained] <- S[constrained]*p[constrained]
        }
        return(-as.vector(S))
    }

   Information <- function(p,gamma,pr) {
        if (optim$constrain) {
            p[constrained] <- exp(p[constrained])
        }
        myp <- lapply(ParPos,function(x) p[x])
        D <- lapply(1:length(myp), function(j) {
            ## K <- ffz[,j]/sffz
            val <- score(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE,model=MODEL)
            apply(val,2,function(x) x*gamma[,j])
        })
        D0 <- matrix(0,nrow(data),length(p))
        for (j in 1:k) D0[,ParPos[[j]]] <- D0[,ParPos[[j]]]+D[[j]]
        if (optim$constrain) {
            for (j in constrained)
                D0[,j] <- D0[,j]*p[j]
        }
        S <- colSums(D0)
        structure(crossprod(D0), grad=S)
    }

    EMstep <- function(p,all=FALSE) {
        thetacur <- p[seq(Npar)]
        gamma <- PosteriorProb(p,constrain=optim$constrain)
        probcur <- colMeans(gamma)
        ## I <- function(p) {
        ##     I <- Information(p,gamma,probcur)
        ##     D <- attr(I, "grad")
        ##     res <- -Inverse(I)
        ##     res <- -I
        ##     attributes(res)$grad <- D
        ##     res
        ## }
        D <- function(p) GradEstep(p,gamma,probcur)
        ## if (optim$newton>0) {
        ##     newpar <- NR(thetacur, gradient=D)
        ## }
        ## if (mean(newpar$gradient^2)>optim$tol) {
        newpar <- nlminb(thetacur,function(p) ObjEstep(p,gamma,probcur), D)
        ## }
        thetacur <- newpar$par
        thetacur0 <- thetacur
        if (optim$constrain) {
            thetacur0[constrained] <- exp(thetacur[constrained])
        }
        p <- c(thetacur,probcur)
        if (all) {
            res <- list(p=p,gamma=gamma,
                        theta=rbind(thetacur0),
                        prob=rbind(probcur))
            return(res)
        }
        return(p)
    }
    em.idx <- match(c("K","method","square","step.min0","step.max0","mstep",
                      "objfn.inc","kr",
                      ##"keep.objfval","convtype",
                      "maxiter","tol","trace"),base::names(optim))
    em.control <- optim[na.omit(em.idx)]
    if (!is.null(em.control$trace)) em.control$trace <- em.control$trace>0

    p <- c(thetacur,probcur)
    opt <- SQUAREM::squarem(p,fixptfn=EMstep,##objfn=negLogLik,
                            control=em.control)
    ## opt2 <- nlminb(opt$par, function(p) negLogLik(p=p), control=list(trace=1))
    val <- EMstep(opt$par,all=TRUE)
    delta <- 1e-6
    if (any(val$prob<delta)) {
        val$prob[val$prob<delta] <- delta
        val$prob <- val$prob/sum(val$prob)
    }

    val <- c(val, list(member=apply(val$gamma,1,which.max),
                       k=k,
                       data=data,
                       parpos=ParPos,
                       multigroup=mg,
                       model=mg$lvm,
                       parname=parname,
                       logLik=NA,
                       opt=opt))
    class(val) <- "lvm.mixture"

    np <- length(coef(val))
    if (is.null(vcov)) {
        val$vcov <- matrix(NA,np,np)
    } else {
        I <- suppressWarnings(information.lvm.mixture(val,type=vcov))
        val$vcov <- tryCatch(Inverse(I),error=function(...) matrix(NA,np,np))
    }
    return(val)
}

###}}} mixture


###{{{ logLik, score, information

##' @export
score.lvm.mixture <- function(x,p=coef(x,full=TRUE),prob,indiv=FALSE,model="normal",...) {
    myp <- lapply(x$parpos,function(x) p[x])
    if (missing(prob)) {
        prob <- tail(p,x$k-1)
        p <- p[seq_len(length(p)-x$k+1)]
    }
    if (length(prob)<x$k)
        prob <- c(prob,1-sum(prob))
    logff <- sapply(seq(x$k), function(j) (logLik(x$multigroup$lvm[[j]],p=myp[[j]],data=x$data,indiv=TRUE,model=model)))
    logplogff <- t(apply(logff,1, function(y) y+log(prob)))
    zmax <- apply(logplogff,1,max)
    logsumpff <- log(rowSums(exp(logplogff-zmax)))+zmax
    aji <- apply(logplogff,2,function(x) exp(x-logsumpff))

    scoref <- lapply(score(x$multigroup,p=p,indiv=TRUE,model=model),
                     function(x) { x[which(is.na(x))] <- 0; x })

    Stheta <- matrix(0,ncol=ncol(scoref[[1]]),nrow=nrow(scoref[[1]]))
    Spi <- matrix(0,ncol=x$k-1,nrow=nrow(Stheta))
    for (j in 1:x$k) {
        Stheta <- Stheta + apply(scoref[[j]],2,function(x) x*aji[,j])
        if (j<x$k)
            Spi[,j] <- aji[,j]/prob[j] - aji[,x$k]/prob[x$k]
    }
    S <- cbind(Stheta,Spi)
    if (!indiv)
        return(colSums(S))
    return(S)
}

##' @export
information.lvm.mixture <- function(x,p=coef(x,full=TRUE),...,type="observed") {
    if (tolower(type)%in%c("obs","observed","hessian")) {
        S <- function(p=p,...) score.lvm.mixture(x,p=p,indiv=FALSE,...)
        I <- -numDeriv::jacobian(S,p)
        return(I)
    }
    S <- score.lvm.mixture(x,p=p,indiv=TRUE,...)
    res <- t(S)%*%S
    attributes(res)$grad <- colSums(S)
    return(res)
}

##' @export
logLik.lvm.mixture <- function(object,p=coef(object,full=TRUE),prob,model="normal",...) {
    if (is.null(object$parpos)) {
        object$parpos <- modelPar(object$multigroup,seq_along(p))$p
    }
    myp <- lapply(object$parpos, function(x) p[x])
    if (missing(prob))
        prob <- tail(p,object$k-1)
    if (length(prob)<object$k)
        prob <- c(prob,1-sum(prob))
    logff <- sapply(seq(object$k), function(j) (logLik(object$multigroup$lvm[[j]],p=myp[[j]],data=object$data,indiv=TRUE,model=model)))
    logplogff <- t(apply(logff,1, function(y) y+log(prob)))
    ## Log-sum-exp (see e.g. NR)
    zmax <- apply(logplogff,1,max)
    logsumpff <- log(rowSums(exp(logplogff-zmax)))+zmax
    loglik <- sum(logsumpff)
    npar <- length(p)
    nobs <- nrow(object$data)
    attr(loglik, "nall") <- nobs
    attr(loglik, "nobs") <- nobs
    attr(loglik, "df") <- npar
    class(loglik) <- "logLik"
    return(loglik)
}

###}}} logLik, score, information

###{{{ vcov

##' @export
vcov.lvm.mixture <- function(object,...) {
    return(object$vcov)
}

###}}}

###{{{ summary/print

##' @export
summary.lvm.mixture <- function(object,type=0,labels=0,...) {
    ppos <- parpos(object$multigroup)
    parname <- attr(ppos,"name")
    naparname <- which(is.na(parname))
    parname[naparname]  <- object$multigroup$name[naparname]
    ParPos <- modelPar(object$multigroup,seq_along(coef(object)))$p

    mm <- object$multigroup$lvm
    p <- coef(object,list=TRUE)
    p0 <- coef(object,prob=FALSE)
    myp <- modelPar(object$multigroup,1:length(p0))$p
    coefs <- list()
    ncluster <- c()
    Coefs <- matrix(NA,ncol=4,nrow=length(parname))
    colnames(Coefs) <- c("Estimate","Std. Error", "Z value", "Pr(>|z|)")
    Types <- rep("other",length(parname))
    Variable  <- rep(NA,length(parname))
    From  <- rep(NA,length(parname))
    Latent <- c()
    for (i in 1:length(mm)) {
        cc.idx <- order(coef(mm[[i]],p=seq_along(p[[i]]),type=2)[,1])

        cc <- coef(mm[[i]],p=p[[i]],vcov=vcov(object)[myp[[i]],myp[[i]]],data=NULL,labels=labels,type=2)
        Latent <- union(Latent,attr(cc,"latent"))
        Coefs[ParPos[[i]],] <- cc[cc.idx,,drop=FALSE]
        Types[ParPos[[i]]] <- attr(cc,"type")[cc.idx]
        Variable[ParPos[[i]]] <- attr(cc, "var")[cc.idx]
        From[ParPos[[i]]] <- attr(cc, "from")[cc.idx]
        cc <- CoefMat(mm[[i]],p=p[[i]],vcov=vcov(object)[myp[[i]],myp[[i]]],data=NULL,labels=labels)
        coefs <- c(coefs, list(cc))
        ncluster <- c(ncluster,sum(object$member==i))
    }
    rownames(Coefs) <- parname

    res <- list(coef=coefs, coefmat=Coefs, coeftype=Types,
                type=type, var=Variable, from=From, latent=Latent,
                ncluster=ncluster, prob=tail(object$prob,1),
                AIC=AIC(object), s2=sum(score(object)^2))
    class(res) <- "summary.lvm.mixture"
    return(res)
}

##' @export
print.summary.lvm.mixture <- function(x,...) {
    if (x$type>0) {
        cc <- x$coefmat
        attr(cc,"type") <- x$coeftype
        attr(cc,"latent") <- x$latent
        attr(cc,"var") <- x$var
        attr(cc,"from") <- x$from
        cat("Mixing parameters:\n")
        cat("  ", paste(as.vector(formatC(x$prob))),"\n")
        print(CoefMat(cc), quote=FALSE)
        return(invisible())
    }
    for (i in 1:length(x$coef)) {
        cat("Cluster ",i," (n=",x$ncluster[i],", Prior=", formatC(x$prob[i]),"):\n",sep="")
        cat(rep("-",50),"\n",sep="")
        print(x$coef[[i]], quote=FALSE)
        if (i<length(x$coef)) cat("\n")
    }
    cat(rep("-",50),"\n",sep="")
    cat("AIC=",x$AIC,"\n")
    cat("||score||^2=",x$s2,"\n")
    invisible()
}

##' @export
print.lvm.mixture <- function(x,...) {
    ## for (i in 1:x$k) {
    ##     cat("Cluster ",i," (pr=",x$prob[i],"):\n",sep="")
    ##     cat(rep("-",50),"\n",sep="")
    ##     print(coef(x,prob=FALSE)[x$parpos[[i]]], quote=FALSE)
    ##     cat("\n")
    ## }
    print(summary(x,type=1,...))
}

###}}}

###{{{ plot

##' @export
plot.lvm.mixture <- function(x,type="l",...) {
    matplot(x$theta,type=type,...)
}

###}}} plot

###{{{ coef

##' @export
coef.lvm.mixture <- function(object,iter,list=FALSE,full=TRUE,prob=FALSE,class=FALSE,label=TRUE,...) {
    N <- nrow(object$theta)
    res <- object$theta
    ## nn <- attr(parpos(object$multigroup),"name")
    ## if (length(ii <- which(is.na(nn)))>0 && label) {
    ##     nn[ii] <- paste0("p",seq_along(ii))
    ## }
    nn <- object$parname
    colnames(res) <- nn
    if (class) return(object$gammas)
    if (list) {
        res <- list()
        for (i in 1:object$k) {
            nn <- coef(object$multigroup$lvm[[i]])
            cc <- coef(object)[object$parpos[[i]]]
            base::names(cc) <- nn
            res <- c(res, list(cc))
        }
        return(res)
    }
    if (full) {
        pp <- object$prob[,seq(ncol(object$prob)-1),drop=FALSE]
        colnames(pp) <- paste0("pr",seq(ncol(pp)))
        res <- cbind(res,pp)
    }
    if (prob) {
        res <- object$prob
    }
    if (missing(iter))
        return(res[N,,drop=TRUE])
    else
        return(res[iter,])
}

###}}} coef

##' @export
model.frame.lvm.mixture <- function(formula,...) {
    return(formula$data)
}

##' @export
iid.lvm.mixture <- function(x,...) {
    bread <- vcov(x)
    structure(t(bread%*%t(score(x,indiv=TRUE))),bread=bread)
}

##' @export
manifest.lvm.mixture <- function(x,...) {
    manifest(x$multigroup,...)
}
