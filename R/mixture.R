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
#' # @param type Type of EM algorithm (standard, classification, stochastic)
#' @param vcov of asymptotic covariance matrix (NULL to omit)
#' @param ... Additional arguments parsed to lower-level functions
#' @author Klaus K. Holst
#' @seealso \code{mvnmix}
#' @keywords models, regression
#' @examples
#'
#' \donttest{
#' m0 <- lvm(list(y~x+z,x~z))
#' distribution(m0,~z) <- binomial.lvm()
#' d <- sim(m0,2000,p=c("y<-z"=2,"y<-x"=1),seed=1)
#'
#' ## unmeasured confounder example
#' m <- baptize(lvm(y~x));
#' covariance(m,~x) <- "v"
#' intercept(m,~x+y) <- NA
#'
#' set.seed(42)
#' M <- mixture(m,k=2,data=d,control=list(trace=1,tol=1e-6))
#' summary(M)
#' lm(y~x,d)
#' estimate(M,"y~x")
#' ## True slope := 1
#' }
#'
#' @export mixture
mixture <- function(x, data, k=length(x),
                    control=list(),
                    vcov="observed",
            ...) {
    ##type=c("standard","CEM","SEM"),     
    MODEL <- "normal"
    ## type <- tolower(type[1])
    ## if (type[1]!="standard") {
    ##     return(mixture0(x,data=data,k=k,control=control,type=type,...))
    ## }

    optim <- list(
        rerun=TRUE,
        ## EM (accelerated):
        K=2, ## K=1, first order approx, slower but some times more stable
        square=TRUE,
        step.max0=1,
        step.min0=1,
        mstep=4,
        kr=1,
        objfn.inc=2,

        keep.objfval=TRUE,
        convtype= "parameter",
        ##convtype = "objfn",
        maxiter=200,
        tol=1e-5,
        trace=1,

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
        optim[names(control)] <- control
    if ("iter.max"%in%names(optim)) optim$maxiter <- optim$iter.max

    if (k==1) {
        if (is.list(x))
            res <- estimate(x[[1]],data,...)
        else
            res <- estimate(x,data,...)
        return(res)
    }
    if (class(x)[1]=="lvm") {
        index(x) <- reindex(x,zeroones=TRUE,deriv=TRUE)
        x <- rep(list(x),k)
    }

    mg <- multigroup(x,rep(list(data),k),fix=FALSE)
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

    if (is.null(optim$start)) {
        constrLogLikS <- function(p) {
            if (optim$constrain) {
                p[constrained] <- exp(p[constrained])
            }
            -logLik(mymodel,c(p=p,rep(1/k,k-1)))
        }

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
        start[constrained] <- exp(start[constrained])
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
    probs <- rbind(probcur);
    thetas <- rbind(thetacur)

    if (optim$constrain) {
        thetacur[constrained] <- log(thetacur[constrained])
    }


    ##  gammas <- list()

    PosteriorProb <- function(pp,priorprob,constrain=FALSE) {
        if (!is.list(pp)) {
            if (constrain) {
                pp[constrained] <- exp(pp[constrained])
            }
            if (missing(priorprob)) priorprob <- pp[seq(Npar+1,length(pp))]
            pp <- lapply(ParPos,function(x) pp[x])
        }
        k <- length(pp)
        logff <- sapply(seq(k), function(j) (logLik(mg$lvm[[j]],p=pp[[j]],data=data,indiv=TRUE,model=MODEL)))
        logplogff <- t(apply(logff,1, function(z) z+log(priorprob)))
        ## Log-sum-exp (see e.g. NR)
        zmax <- apply(logplogff,1,max)
        logsumpff <- log(rowSums(exp(logplogff-zmax)))+zmax
        gamma <- exp(apply(logplogff,2,function(y) y - logsumpff)) ## Posterior class probabilities
        return(gamma)
    }

    myObj <- function(p,gamma) {
        if (optim$constrain) {
            p[constrained] <- exp(p[constrained])
        }
        if (missing(gamma)) {
            gamma <- PosteriorProb(p)
        }
        myp <- lapply(ParPos,function(x) p[x])
        prob <- p[seq(Npar+1,length(p))]

        logff <- sapply(1:length(myp), function(j) (logLik(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE,model=MODEL)))
        logplogff <- t(apply(logff,1, function(y) y+log(prob)))
        zmax <- apply(logplogff,1,max)
        logsumpff <- log(rowSums(exp(logplogff-zmax)))+zmax
        loglik <- sum(logsumpff)
        return(-loglik)
    }

    Scoring <- function(p,gamma) {
        p.orig <- p
        if (optim$constrain) {
            p[constrained] <- exp(p[constrained])
        }
        probcur <- colMeans(gamma)
        myp <- lapply(ParPos,function(x) p[x])
        D <- lapply(seq_along(myp), function(j) gamma[,j]*score(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE,model=MODEL))
        D0 <- matrix(0,nrow(data),length(p))
        for (j in 1:k) D0[,ParPos[[j]]] <- D0[,ParPos[[j]]]+D[[j]]
        S <- colSums(D0)
        if (optim$constrain) {
            S[constrained] <- S[constrained]*p[constrained]
        }
        I <- lapply(1:k, function(j)
            probcur[j]*information(mg$lvm[[j]],p=myp[[j]],n=nrow(data),data=data,model=MODEL))
        I0 <- matrix(0,length(p),length(p))
        for (j in 1:k) {
            I0[ParPos[[j]],ParPos[[j]]] <- I0[ParPos[[j]],ParPos[[j]]] + I[[j]]
        }
        if (optim$constrain) {
            I0[constrained,-constrained] <- apply(I0[constrained,-constrained,drop=FALSE],2,function(x) x*p[constrained]);
            I0[-constrained,constrained] <- t(I0[constrained,-constrained])
            if (length(constrained)==1)
                I0[constrained,constrained] <- I0[constrained,constrained]*p[constrained]^2 + S[constrained]
            else
                I0[constrained,constrained] <- I0[constrained,constrained]*outer(p[constrained],p[constrained]) + diag(S[constrained])
        }
        ##    print(paste(S,collapse=","))
        if (optim$stabil) {
            ##      I0 <- I0+S%*%t(S)
            if (optim$lambda>0)
                sigma <- optim$lambda
            else
                sigma <- (t(S)%*%S)[1]^0.5
            I0 <- I0+optim$gamma2*(sigma)*diag(nrow(I0))
        }
        p.orig + as.vector(optim$gamma*Inverse(I0)%*%S)
        ##   p.orig + Inverse(I0+optim$lambda*diag(nrow(I0)))%*%S
    }
    ## env <- new.env()
    ## assign("mg",mg,env)
    ## assign("k",k,env)
    ## assign("data",data,env)
    ## assign("MODEL",MODEL,env)
    ## assign("optim",optim,env)
    ## assign("parpos",parpos,env)
    ## assign("constrained",constrained,env)


    EMstep <- function(p,all=FALSE) {
        thetacur0 <- thetacur <- p[seq(Npar)]
        gamma <- PosteriorProb(p,constrain=optim$constrain)
        probcur <- colMeans(gamma)
        count2 <- 0
        for (jj in 1:optim$newton) {
            count2 <- count2+1
            oldpar <- thetacur
            thetacur <- Scoring(thetacur,gamma)
            if (frobnorm(oldpar-thetacur)<optim$delta) break;
        }
        theteacur0 <- thetacur
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

    ## A few steps with the standard EM
    ## p0 <- p
    ## for (i in 1:5) {
    ##     p0 <- EMstep(p0)
    ## }

    ## sqem.idx <- match(c("K","method","square","step.min0","step.max0","mstep",
    ##                         "objfn.inc","kr","tol","maxiter","trace"),names(optim))
    ## em.control <- optim[na.omit(sqem.idx)]
    ## run.idx <- match(c("keep.objfval","convtype","maxiter","tol","trace"),names(optim))
    ## control.run <- optim[na.omit(run.idx)]

    em.idx <- match(c("K","method","square","step.min0","step.max0","mstep",
                      "objfn.inc","kr",
                      ##"keep.objfval","convtype",
                      "maxiter","tol","trace"),names(optim))
    em.control <- optim[na.omit(em.idx)]
    if (!is.null(em.control$trace)) em.control$trace <- em.control$trace>0

    p <- c(thetacur,probcur)
    ## opt <- turboEM::turboem(p,fixptfn=EMstep,
    ##                         objfn=myObj,
    ##                         method="squarem",
    ##                         ##method = c("em","squarem","pem","decme","qn")
    ##                         control.method=list(em.control),
    ##                         control.run=control.run)
    opt <- SQUAREM::squarem(p,fixptfn=EMstep,#objfn=myObj,
                            control=em.control)
    ## opt <- SQUAREM::fpiter(p,fixptfn=EMstep,control=list(maxiter=100,trace=1))

    val <- EMstep(opt$par,all=TRUE)
    delta <- 1e-4
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
                       logLik=NA,
                       opt=opt))
    class(val) <- "lvm.mixture"

    ## if (optim$rerun) {
    ##     p0 <- coef(val)
    ##     f <- function(p) -logLik(val,p)
    ##     g <- function(p) -score(val,p)
    ##     suppressWarnings(opt <- tryCatch(ucminf(p0,f,g),error=function(x) opt))
    ## }
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
  ## if (missing(prob))
  ##   prob <- coef(x,prob=TRUE)
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
summary.lvm.mixture <- function(object,labels=0,...) {
  mm <- object$multigroup$lvm
  p <- coef(object,list=TRUE)
  p0 <- coef(object,prob=FALSE)
  myp <- modelPar(object$multigroup,1:length(p0))$p
  coefs <- list()
  ncluster <- c()
  for (i in 1:length(mm)) {
    cc <- CoefMat(mm[[i]],p=p[[i]],vcov=vcov(object)[myp[[i]],myp[[i]]],data=NULL,labels=labels)
    coefs <- c(coefs, list(cc))
    ncluster <- c(ncluster,sum(object$member==i))
  }
  res <- list(coef=coefs,ncluster=ncluster,prob=tail(object$prob,1),
              AIC=AIC(object),s2=sum(score(object)^2))
  class(res) <- "summary.lvm.mixture"
  return(res)
}

##' @export
print.summary.lvm.mixture <- function(x,...) {
  space <- paste(rep(" ",12),collapse="")
  for (i in 1:length(x$coef)) {
    cat("Cluster ",i," (n=",x$ncluster[i],", Prior=", formatC(x$prob[i]),"):\n",sep="")
    cat(rep("-",50),"\n",sep="")
    print(x$coef[[i]], quote=FALSE)
    if (i<length(x$coef)) cat("\n")
  }
  cat(rep("-",50),"\n",sep="")
  cat("AIC=",x$AIC,"\n")
  cat("||score||^2=",x$s2,"\n")
  invisible(par)  
}

##' @export
print.lvm.mixture <- function(x,...) {
  space <- paste(rep(" ",12),collapse="")
  for (i in 1:x$k) {
    cat("Cluster ",i," (n=",sum(x$member==i),"):\n",sep="")
    cat(rep("-",50),"\n",sep="")
    print(coef(x,prob=FALSE)[x$parpos[[i]]], quote=FALSE)
    cat("\n")   
  }
  invisible(par)
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
  nn <- attr(parpos(object$multigroup),"name")
  if (length(ii <- which(is.na(nn)))>0 && label) {
      nn[ii] <- paste0("p",seq_along(ii))
  }
  colnames(res) <- nn
  if (class) return(object$gammas)
  if (list) {      
      res <- list()
      for (i in 1:object$k) {
          nn <- coef(object$multigroup$lvm[[i]])
          cc <- coef(object)[object$parpos[[i]]]
          names(cc) <- nn
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
