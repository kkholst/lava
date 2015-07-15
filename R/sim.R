##' Simulate model
##'
##' Simulate data from a general SEM model including non-linear effects and
##' general link and distribution of variables.
##'
##' @aliases sim sim.lvmfit sim.lvm
##' simulate.lvmfit simulate.lvm
##' transform<- transform<-.lvm transform.lvm
##' functional functional<-  functional.lvm functional<-.lvm
##' distribution distribution distribution<- distribution.lvm distribution<-.lvm
##' heavytail heavytail<-
##' weibull.lvm
##' binomial.lvm
##' poisson.lvm
##' uniform.lvm
##' normal.lvm
##' lognormal.lvm
##' gaussian.lvm
##' probit.lvm
##' logit.lvm
##' pareto.lvm
##' student.lvm
##' chisq.lvm
##' coxGompertz.lvm
##' coxWeibull.lvm
##' coxExponential.lvm
##' aalenExponential.lvm
##' Gamma.lvm gamma.lvm
##' loggamma.lvm
##' categorical categorical<-
##' threshold.lvm
##' ones.lvm
##' sequence.lvm
##' @usage
##' \method{sim}{lvm}(x, n = 100, p = NULL, normal = FALSE, cond = FALSE,
##' sigma = 1, rho = 0.5, X, unlink=FALSE, ...)
##' @param x Model object
##' @param n Number of simulated values/individuals
##' @param p Parameter value (optional)
##' @param normal Logical indicating whether to simulate data from a
##' multivariate normal distribution conditional on exogenous variables hence
##' ignoring functional/distribution definition
##' @param cond for internal use
##' @param sigma Default residual variance (1)
##' @param rho Default covariance parameter (0.5)
##' @param X Optional matrix of covariates
##' @param unlink Return Inverse link transformed data
##' @param \dots Additional arguments to be passed to the low level functions
##' @author Klaus K. Holst
##' @keywords models datagen regression
##' @export
##' @examples
##' ##################################################
##' ## Logistic regression
##' ##################################################
##' m <- lvm(y~x+z)
##' regression(m) <- x~z
##' distribution(m,~y+z) <- binomial.lvm("logit")
##' d <- sim(m,1e3)
##' head(d)
##'
##' e <- estimate(m,d,estimator="glm")
##' e
##' ## Simulate a few observation from estimated model
##' sim(e,n=5)
##'
##' ##################################################
##' ## Poisson
##' ##################################################
##' distribution(m,~y) <- poisson.lvm()
##' d <- sim(m,1e4,p=c(y=-1,"y~x"=2,z=1))
##' head(d)
##' estimate(m,d,estimator="glm")
##' mean(d$z); lava:::expit(1)
##'
##' summary(lm(y~x,sim(lvm(y[1:2]~4*x),1e3)))
##'
##' ##################################################
##' ### Gamma distribution
##' ##################################################
##' m <- lvm(y~x)
##' distribution(m,~y+x) <- list(Gamma.lvm(shape=2),binomial.lvm())
##' intercept(m,~y) <- 0.5
##' d <- sim(m,1e4)
##' summary(g <- glm(y~x,family=Gamma(),data=d))
##' \dontrun{MASS::gamma.shape(g)}
##'
##' args(lava::Gamma.lvm)
##' distribution(m,~y) <- Gamma.lvm(shape=2,log=TRUE)
##' sim(m,10,p=c(y=0.5))[,"y"]
##'
##' ##################################################
##' ### Transform
##' ##################################################
##'
##' m <- lvm()
##' transform(m,xz~x+z) <- function(x) x[1]*(x[2]>0)
##' regression(m) <- y~x+z+xz
##' d <- sim(m,1e3)
##' summary(lm(y~x+z + x*I(z>0),d))
##'
##'
##' ##################################################
##' ### Non-random variables
##' ##################################################
##' m <- lvm()
##' distribution(m,~x+z+v+w) <- list(sequence.lvm(0,5),## Seq. 0 to 5 by 1/n
##'                                ones.lvm(),       ## Vector of ones
##'                                ones.lvm(0.5),    ##  0.8n 0, 0.2n 1
##'                                ones.lvm(interval=list(c(0.3,0.5),c(0.8,1))))
##' sim(m,10)
##'
##'
##' ##################################################
##' ### Cox model
##' ### piecewise constant hazard
##' ################################################
##'
##' m <- lvm(t~x)
##' rates <- c(1,0.5); cuts <- c(0,5)
##' ## Constant rate: 1 in [0,5), 0.5 in [5,Inf)
##' distribution(m,~t) <- coxExponential.lvm(rate=rates,timecut=cuts)
##'
##'
##' \dontrun{
##'     d <- sim(m,2e4,p=c("t~x"=0.1)); d$status <- TRUE
##'     plot(timereg::aalen(survival::Surv(t,status)~x,data=d,
##'                         resample.iid=0,robust=0),spec=1)
##'     L <- approxfun(c(cuts,max(d$t)),f=1,
##'                    cumsum(c(0,rates*diff(c(cuts,max(d$t))))),
##'                    method="linear")
##'     curve(L,0,100,add=TRUE,col="blue")
##' }
##'
##'
##' ##################################################
##' ### Cox model
##' ### piecewise constant hazard, gamma frailty
##' ##################################################
##'
##' m <- lvm(y~x+z)
##' rates <- c(0.3,0.5); cuts <- c(0,5)
##' distribution(m,~y+z) <- list(coxExponential.lvm(rate=rates,timecut=cuts),
##'                              loggamma.lvm(rate=1,shape=1))
##' \dontrun{
##'     d <- sim(m,2e4,p=c("y~x"=0,"y~z"=0)); d$status <- TRUE
##'     plot(timereg::aalen(survival::Surv(y,status)~x,data=d,
##'                         resample.iid=0,robust=0),spec=1)
##'     L <- approxfun(c(cuts,max(d$y)),f=1,
##'                    cumsum(c(0,rates*diff(c(cuts,max(d$y))))),
##'                    method="linear")
##'     curve(L,0,100,add=TRUE,col="blue")
##' }
##'
##' ## Equivalent via transform (here with Aalens additive hazard model)
##' m <- lvm(y~x)
##' distribution(m,~y) <- aalenExponential.lvm(rate=rates,timecut=cuts)
##' distribution(m,~z) <- Gamma.lvm(rate=1,shape=1)
##' transform(m,t~y+z) <- prod
##' sim(m,10)
##'
##' ## Shared frailty
##' m <- lvm(c(t1,t2)~x+z)
##' rates <- c(1,0.5); cuts <- c(0,5)
##' distribution(m,~y) <- aalenExponential.lvm(rate=rates,timecut=cuts)
##' distribution(m,~z) <- loggamma.lvm(rate=1,shape=1)
##' \dontrun{
##' mets::fast.reshape(sim(m,100),varying="t")
##' }
##'
##' ##################################################
##' ### General multivariate distributions
##' ##################################################
##'
##' \dontrun{
##' m <- lvm()
##' distribution(m,~y1+y2,oratio=4) <- VGAM::rbiplackcop
##' ksmooth2(sim(m,1e4),rgl=FALSE,theta=-20,phi=25)
##'
##' m <- lvm()
##' distribution(m,~z1+z2,"or1") <- VGAM::rbiplackcop
##' distribution(m,~y1+y2,"or2") <- VGAM::rbiplackcop
##' sim(m,10,p=c(or1=0.1,or2=4))
##' }
##'
##' m <- lvm()
##' distribution(m,~y1+y2+y3,TRUE) <- function(n,...) rmvn(n,sigma=diag(3)+1)
##' var(sim(m,100))
##'
##' ## Syntax also useful for univariate generators, e.g.
##' m <- lvm(y~x+z)
##' distribution(m,~y,TRUE) <- function(n) rnorm(n,mean=1000)
##' sim(m,5)
##' distribution(m,~y,"m1",0) <- rnorm
##' sim(m,5)
##' sim(m,5,p=c(m1=100))
##'
##' ##################################################
##' ### Categorical predictor
##' ##################################################
##'
##' m <- lvm()
##' ## categorical(m,K=3) <- "v"
##' categorical(m,labels=c("A","B","C")) <- "v"
##'
##' regression(m,additive=FALSE) <- y~v
##' \dontrun{
##' plot(y~v,sim(m,1000,p=c("y~v:2"=3)))
##' }
##'
##' m <- lvm()
##' categorical(m,labels=c("A","B","C"),p=c(0.5,0.3)) <- "v"
##' regression(m,additive=FALSE,beta=c(0,2,-1)) <- y~v
##' ## equivalent to:
##' ## regression(m,y~v,additive=FALSE) <- c(0,2,-1)
##' regression(m,additive=FALSE,beta=c(0,4,-1)) <- z~v
##' table(sim(m,1e4)$v)
##' glm(y~v, data=sim(m,1e4))
##' glm(y~v, data=sim(m,1e4,p=c("y~v:1"=3)))
"sim" <- function(x,...) UseMethod("sim")

##' @export
sim.lvmfit <- function(x,n=nrow(model.frame(x)),p=pars(x),xfix=TRUE,...) {
    m <- Model(x)
    if ((nrow(model.frame(x))==n) & xfix) {
        X <- exogenous(x)
        mydata <- model.frame(x)
        for (pred in X) {
            distribution(m, pred) <- list(mydata[,pred])
        }
    }
    sim(m,n=n,p=p,...)
}

##' @export
sim.lvm <- function(x,n=100,p=NULL,normal=FALSE,cond=FALSE,sigma=1,rho=.5,
                    X,unlink=FALSE,...) {
    if (!missing(X)) {
        n <- nrow(X)
    }
    xx <- exogenous(x)
    if (!is.null(p)) {
        i1 <- na.omit(c(match(names(p),xx),
                        match(names(p),paste0(xx,lava.options()$symbol[2],xx))))
        if (length(i1)>0) covariance(x) <- xx[i1]
    }
    ##  index(x) <- reindex(x)
    vv <- vars(x)
    nn <- setdiff(vv,parameter(x))
    mu <- unlist(lapply(x$mean, function(l) ifelse(is.na(l)|is.character(l),0,l)))
    xf <- intersect(unique(parlabels(x)),xx)
    xfix <- c(randomslope(x),xf); if (length(xfix)>0) normal <- FALSE

    if (length(p)!=(index(x)$npar+index(x)$npar.mean+index(x)$npar.ex) | is.null(names(p))) {
        nullp <- is.null(p)
        p0 <- p
        ep <- NULL
        ei <- which(index(x)$e1==1)
        if (length(ei)>0)
            ep <- unlist(x$expar)[ei]
        p <- c(rep(1, index(x)$npar+index(x)$npar.mean),ep)
        p[seq_len(index(x)$npar.mean)] <- 0
        p[index(x)$npar.mean + variances(x)] <- sigma
        p[index(x)$npar.mean + offdiags(x)] <- rho
        if (!nullp) {
            c1 <- coef(x,mean=TRUE,fix=FALSE)
            c2 <- coef(x,mean=TRUE,fix=FALSE,labels=TRUE)
            idx1 <- na.omit(match(names(p0),c1))
            idx11 <- na.omit(match(names(p0),c2))
            idx2 <- na.omit(which(names(p0)%in%c1))
            idx22 <- na.omit(which(names(p0)%in%c2))
            if (length(idx1)>0 && !is.na(idx1))
                p[idx1] <- p0[idx2]
            if (length(idx11)>0 && !is.na(idx11))
                p[idx11] <- p0[idx22]
        }
    }
    M <- modelVar(x,p,data=NULL)
    A <- M$A; P <- M$P
    if (!is.null(M$v)) mu <- M$v

    ## dontsim <- names(distribution(x))[unlist(lapply(distribution(x),function(x) identical(x,NA)))]
    PP <- with(svd(P), v%*%diag(sqrt(d),nrow=length(d))%*%t(u))
    mdist <- distribution(x,multivariate=TRUE)$var
    mdistnam <- names(mdist)
    mii <- match(mdistnam,vars(x))

    if (length(distribution(x))>0 ) {
        ii <- match(names(distribution(x)),vars(x))
        jj <- setdiff(seq(ncol(P)),c(ii,mii))
        E <- matrix(0,ncol=ncol(P),nrow=n)
        if (length(jj)>0)
            E[,jj] <-  matrix(rnorm(length(jj)*n),ncol=length(jj))%*%PP[jj,jj,drop=FALSE]
    } else {
        E <- matrix(rnorm(ncol(P)*n),ncol=ncol(P))%*%PP  ## Error term for conditional normal distributed variables
    }

    if (length(mdistnam)>0) {
        fun <- distribution(x,multivariate=TRUE)$fun
        for (i in seq_along(fun)) {
            mv <- names(which(unlist(mdist)==i))
            ii <- match(mv,vars(x))
            E[,ii] <- distribution(x,multivariate=TRUE)$fun[[i]](n,p=p,object=x,...)
        }
    }

    ## Simulate exogenous variables (covariates)
    res <- matrix(0,ncol=length(nn),nrow=n); colnames(res) <- nn

    vartrans <- names(x$attributes$transform)
    xx <- unique(c(exogenous(x, latent=FALSE, index=TRUE),xfix))
    xx <- setdiff(xx,vartrans)

    X.idx <- match(xx,vv)
    res[,X.idx] <- t(mu[X.idx]+t(E[,X.idx]))
    if (missing(X)) {
        if (!is.null(xx) && length(xx)>0)
            for (i in seq_along(xx)) {
                mu.x <- mu[X.idx[i]]
                dist.x <- distribution(x,xx[i])[[1]]
                if (is.list(dist.x) && is.function(dist.x[[1]])) dist.x <- dist.x[[1]]
                if (is.list(dist.x)) {
                    dist.x <- dist.x[[1]]
                    if (length(dist.x)==1) dist.x <- rep(dist.x,n)
                }
                if (is.function(dist.x)) {
                    res[,X.idx[i]] <- dist.x(n=n,mu=mu.x,var=P[X.idx[i],X.idx[i]])
                } else {
                    if (is.null(dist.x) || is.na(dist.x)) {
                    } else {
                        if (length(dist.x)!=n) stop("'",vv[X.idx[i]], "' fixed at length ", length(dist.x)," != ",n,".")
                        res[,X.idx[i]] <- dist.x ## Deterministic
                    }
                }
            }
    } else {
        res[,X.idx] <- X[,xx]
    }
    simuled <- c(xx)
    resunlink <- NULL
    if (unlink) {
        resunlink <- res
    }

    if ( normal | ( is.null(distribution(x)) & is.null(functional(x)) & is.null(constrain(x))) ) {
        if(cond) { ## Simulate from conditional distribution of Y given X
            mypar <- pars(x,A,P,mu)
            Ey.x <- predict(x, mypar, data.frame(res))
            Vy.x <- attributes(Ey.x)$cond.var
            PP <- with(svd(Vy.x), v%*%diag(sqrt(d),nrow=length(d))%*%t(u))
            yy <- Ey.x + matrix(n*ncol(Vy.x),ncol=ncol(Vy.x))%*%PP
            res <- cbind(yy, res[,xx]); colnames(res) <- c(colnames(Vy.x),xx)
            return(res)
        }
        ## Simulate from sim. distribution (Y,X) (mv-normal)
        I <- diag(length(nn))
        IAi <- Inverse(I-t(A))
        colnames(E) <- vv
        dd <- t(apply(heavytail.sim.hook(x,E),1,function(x) x+mu))
        res <- dd%*%t(IAi)
        colnames(res) <- vv
    } else {


        xconstrain.idx <- unlist(lapply(lapply(constrain(x),function(z) attributes(z)$args),function(z) length(intersect(z,index(x)$manifest))>0))
        xconstrain <- intersect(unlist(lapply(constrain(x),function(z) attributes(z)$args)),index(x)$manifest)

        ##    if (!all(xconstrain %in% index(x)$exogenous)) warning("Non-linear constraint only allowed via covariates")
        if (length(xconstrain)>0)
          for (i in which(xconstrain.idx)) {
            ff <- constrain(x)[[i]]
            myargs <- attributes(ff)$args
            D <- matrix(0,n,length(myargs))
            for (j in seq_len(ncol(D))) {
              if (myargs[j]%in%xconstrain)
                D[,j] <- res[,myargs[j]]
              else
                D[,j] <- M$parval[[myargs[j]]]
            }
            ##res[,names(xconstrain.idx)[i]] <- apply(D,1,ff)
            res <- cbind(res, apply(D,1,ff)); colnames(res)[ncol(res)] <- names(xconstrain.idx)[i]
          }

        xconstrain.par <- names(xconstrain.idx)[xconstrain.idx]
        covparnames <- unique(as.vector(covariance(x)$labels))

        if (any(xconstrain.par%in%covparnames)) {
            mu0 <- rep(0,ncol(P))
            P0 <- P
            E <- t(sapply(seq_len(n),function(idx) {
                for (i in intersect(xconstrain.par,covparnames)) {
                    P0[covariance(x)$labels==i] <- res[idx,i]
                }
                ##        return(rmvnorm(1,mu0,P0))
                PP <- with(svd(P0), v%*%diag(sqrt(d),nrow=length(d))%*%t(u))
                return(mu0+rbind(rnorm(ncol(P0)))%*%PP)
            }))
        } else {
        }
        colnames(E) <- vv
        E <- heavytail.sim.hook(x,E)

        ## Non-linear regression components
        xconstrain <- c()
        for (i in seq_len(length(constrain(x)))) {
            z <- constrain(x)[[i]]
            xx <- intersect(attributes(z)$args,manifest(x))
            if (length(xx)>0) {
                warg <- setdiff(attributes(z)$args,xx)
                wargidx <- which(attributes(z)$args%in%warg)
                exoidx <- which(attributes(z)$args%in%xx)
                parname <- names(constrain(x))[i]
                y <- names(which(unlist(lapply(intercept(x),function(x) x==parname))))
                el <- list(i,y,parname,xx,exoidx,warg,wargidx,z)
                names(el) <- c("idx","endo","parname","exo","exoidx","warg","wargidx","func")
                xconstrain <- c(xconstrain,list(el))
            }
        }

        for (i in intersect(exogenous(x),names(x$constrainY))) {
            cc <- x$constrainY[[i]]
            args <- cc$args
            args <- if (is.null(args) || length(args)==0) res[,i] else res[,args]
            res[,i] <- cc$fun(args,p,...)
        }
        yconstrain <- unlist(lapply(xconstrain,function(x) x$endo))

        res <- data.frame(res)
        if (length(vartrans)>0) {
            parvals <- parpos(x)$parval
            parvalsnam <- setdiff(names(parvals),xx)
            if (length(parvalsnam)>0) {
                Parvals <- p[unlist(parvals)];
                res <- cbind(res,
                             cbind(rep(1,nrow(res)))%x%rbind(Parvals))
                colnames(res)[seq(length(Parvals))+ncol(res)-length(Parvals)] <-
                    names(parvals)
            }
        }

        leftovers <- c()
        while (length(simuled)<length(nn)) {
            leftoversPrev <- leftovers
            leftovers <- setdiff(nn,simuled)

            if (!is.null(leftoversPrev) && length(leftoversPrev)==length(leftovers)) stop("Infinite loop (probably problem with 'transform' call in model: Outcome variable should not affect other variables in the model).")

            for (i in leftovers) {
                if (i%in%vartrans) {
                    xtrans <- x$attributes$transform[[i]]$x
                    if (all(xtrans%in%c(simuled,names(parvals))))  {
                        suppressWarnings(yy <- with(x$attributes$transform[[i]],fun(res[,xtrans])))
                        if (length(yy) != NROW(res)) { ## apply row-wise
                            res[,i] <- with(x$attributes$transform[[i]],apply(res[,xtrans,drop=FALSE],1,fun))
                        } else {
                            colnames(yy) <- NULL
                            res[,i] <- yy
                        }
                        simuled <- c(simuled,i)
                    }
                } else {

                    ipos <- which(i%in%yconstrain)
                    if (length(ipos)==0 || all(xconstrain[[ipos]]$exo%in%simuled)) {
                        pos <- match(i,vv)
                        relations <- colnames(A)[A[,pos]!=0]
                        simvars <- x$attributes$simvar[[i]]

                        if (all(c(relations,simvars)%in%simuled)) { ## Only depending on already simulated variables
                            if (x$mean[[pos]]%in%xconstrain.par) {
                                mu.i <- res[,x$mean[[pos]] ]
                            } else {
                                mu.i <- mu[pos]
                            }
                            if (length(ipos)>0) {
                                pp <- unlist(M$parval[xconstrain[[ipos]]$warg])
                                myidx <- with(xconstrain[[i]],order(c(wargidx,exoidx)))
                                mu.i <- mu.i + with(xconstrain[[ipos]],
                                                    apply(res[,exo,drop=FALSE],1,
                                                          function(x) func(
                                                              unlist(c(pp,x))[myidx])))
                            }

                            for (From in relations) {
                                f <- functional(x,i,From)[[1]]
                                if (!is.function(f))
                                    f <- function(x,...) x
                                reglab <- regfix(x)$labels[From,pos]
                                if (reglab%in%c(xfix,xconstrain.par)) {
                                    if (is.function(f)) {
                                        if (length(formals(f))>1) {
                                            mu.i <- mu.i + res[,reglab]*f(res[,From],p)
                                        } else {
                                            mu.i <- mu.i + res[,reglab]*f(res[,From])
                                        }
                                    } else  mu.i <- mu.i + res[,reglab]*res[,From]
                                }
                                else {
                                    if (is.function(f)) {
                                        if (length(formals(f))>1) {
                                            mu.i <- mu.i + A[From,pos]*f(res[,From],p)
                                        } else {
                                            mu.i <- mu.i + A[From,pos]*f(res[,From])
                                        }
                                    } else  mu.i <- mu.i + A[From,pos]*res[,From]
                                }
                            }
                            dist.i <- distribution(x,i)[[1]]
                            if (!is.function(dist.i)) {
                                res[,pos] <- mu.i + E[,pos]
                                if (unlink)
                                    resunlink[,pos] <- res[,pos]
                            }
                            else {
                                if (length(simvars)>0) { ## Depends on mu and also on other variables (e.g. time-depending effect)
                                    if (length(mu.i)==1) mu.i <- rep(mu.i,n)
                                    mu.i <- cbind("m0"=mu.i,res[,simvars,drop=FALSE])
                                }
                                res[,pos] <- dist.i(n=n,mu=mu.i,var=P[pos,pos])
                                if (unlink)
                                    resunlink[,pos] <- mu.i
                            }
                            if (i%in%names(x$constrainY)) {
                                cc <- x$constrainY[[i]]
                                args <- cc$args
                                args <- if (is.null(args) || length(args)==0) res[,pos] else res[,args]
                                res[,pos] <- cc$fun(args,p,...)
                            }
                            simuled <- c(simuled,i)
                        }
                    }
                }
            }
        }
        res <- res[,nn,drop=FALSE]
    }

    res <- as.data.frame(res)
    myhooks <- gethook("sim.hooks")
    for (f in myhooks) {
        res <- do.call(f, list(x=x,data=res,p=p,modelpar=M))
    }
    if (unlink) res <- resunlink

    res <- as.data.frame(res)
    self <- x$attributes$selftransform
    for (v in names(self)) {
        res[,v] <- self[[v]](res[,v])
    }
    return(res)
}



##' @export
simulate.lvm <- function(object,nsim,seed=NULL,...) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    sim(object,nsim,...)
}

##' @export
simulate.lvmfit <- function(object,nsim,seed=NULL,...) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    sim(object,nsim,...)
}

##' Wrapper function for mclapply
##'
##' @export
##' @param x function or 'sim' object
##' @param R Number of replications
##' @param colnames Optional column names
##' @param messages Messages
##' @param mc.cores Number of cores to use
##' @param blocksize Split computations in blocks
##' @param ... Additional arguments to mclapply
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
##' 
##' val <- sim(onerun,R=10,b0=1,messages=0,mc.cores=1)
##' val
##' val <- sim(val,R=40,b0=1,mc.cores=1) ## append results
##' 
##' summary(val,estimate=c(1,1),confint=c(3,4,6,7),true=c(1,1))
##' summary(val,estimate=c(1,1),se=c(2,5),names=c("Model","Sandwich"))
##' 
##' if (interactive()) {
##'     plot(val[1:20,c(1:4)],vline=c(1,NA,NA,NA),yax.flip=TRUE,start=5)
##' }
sim.default <- function(x=NULL,R=100,f=NULL,colnames=NULL,messages=1L,mc.cores=parallel::detectCores(),blocksize=2L*mc.cores,...) {
    requireNamespace("parallel",quietly=TRUE)
    olddata <- NULL
    if (inherits(x,c("data.frame","matrix"))) olddata <- x
    if (inherits(x,"sim")) {
        x <- attr(x,"f")
        if (!is.null(f)) x <- f
    } else {
        if (!is.null(f)) x <- f
        if (!is.function(x)) stop("Expected a function or 'sim' object.")
    }
    if (is.null(x)) stop("Must give new function argument 'f'.")
    res <- val <- NULL
    mycall <- match.call()
    on.exit({
        if (messages>0) close(pb)
        if (is.null(colnames) && !is.null(val)) {
            if (is.matrix(val[[1]])) {
                colnames <- base::colnames(val[[1]])
            } else {
                colnames <- names(val[[1]])
            }
        }
        base::colnames(res) <- colnames
        if (!is.null(olddata)) res <- rbind(olddata,res)
        attr(res,"call") <- mycall
        attr(res,"f") <- x
        class(res) <- c("sim","matrix")
        if (idx.done<R) {
            res <- res[seq(idx.done),,drop=FALSE]
        }
        return(res)
    })
    nfolds <- max(1,round(R/blocksize))
    idx <- split(1:R,sort((1:R)%%nfolds))
    idx.done <- 0
    count <- 0
    if (messages>0) pb <- txtProgressBar(style=3,width=40)
    for (ii in idx) {
        count <- count+1
        val <- parallel::mclapply(ii,x,mc.cores=mc.cores,...)
        if (messages>0) ##getTxtProgressBar(pb)<(i/R)) {
            setTxtProgressBar(pb, count/length(idx))

        if (is.null(res)) {
            res <- matrix(NA,ncol=length(val[[1]]),nrow=R)
        }
        res[ii,] <- Reduce(rbind,val)
        idx.done <- max(ii)
    }
}

##' @export
"[.sim" <- function (x, i, j, drop = FALSE) {
    atr <- attributes(x)
    class(x) <- "matrix"
    x <- NextMethod("[",drop=FALSE)
    atr.keep <- "call"
    if (missing(j)) atr.keep <- c(atr.keep,"f")
    attributes(x)[atr.keep] <- atr[atr.keep]
    class(x) <- c("sim","matrix")
    x
}

##' @export
print.sim <- function(x,...) {
    attr(x,"f") <- attr(x,"call") <- NULL
    class(x) <- "matrix"
    print(x,...)
}

##' @export
plot.sim <- function(x,vline=NULL,lty.vline=2,
                     start=1,end=nrow(x),
                     idx=seq(ncol(x)),
                     main="Running Average",
                     plot.type=c("multiple","single"),xlab="",...) {
    if (!requireNamespace("zoo",quietly=TRUE)) stop("zoo package required")
    my.panel <- function(x, ..., pf = parent.frame()) {
        lines(x, ...)
        args <- list(...)
        args$type <- NULL
        args$lty <- lty.vline
        args[1] <- NULL
        args$h <- vline[with(pf,panel.number)]
        if (!is.null(vline)) do.call(abline,args)
    }
    if (!is.null(vline) && length(idx)>vline) vline <- rep(vline,ceiling(length(idx)/length(vline)))
    val <- apply(x[,idx,drop=FALSE],2,function(z) cumsum(z)/seq(length(z)))
    val <- stats::window(zoo::zoo(val),start=start,end=end)
    plot(val,xlab=xlab,plot.type=plot.type,panel=my.panel,main=main,...)
}

##' @export
summary.sim <- function(object,estimate=NULL,se=NULL,confint=NULL,true=NULL,fun,names=NULL,...) {
    if (missing(fun)) fun <- function(x) {
        pp <- c(.025,.5,.975)
        res <- c(mean(x,na.rm=TRUE),sd(x,na.rm=TRUE),quantile(x,c(0,pp,1),na.rm=TRUE),
                 mean(is.na(x)))
        names(res) <- c("Mean","SD","Min",paste0(pp*100,"%"),"Max","Missing")
        res
    }
    if (is.null(estimate) && is.null(confint)) return(apply(object,2,fun))

    est <- apply(object[,estimate,drop=FALSE],2,
                 function(x) c(Mean=mean(x,na.rm=TRUE),Missing=mean(is.na(x)),SD=sd(x,na.rm=TRUE)))
    if (!is.null(true)) {
        if (length(true)!=length(estimate)) stop("'true' should be of same length as 'estimate'.")
        est <- rbind(rbind(True=true),rbind(Bias=true-est["Mean",]),
                     rbind(RMSE=((true-est["Mean",])^2+(est["SD",])^2)^.5),
                     est)
    }
    if (!is.null(se)) {
        if (length(se)!=length(estimate)) stop("'se' should be of same length as 'estimate'.")
        est <- rbind(est, SE=apply(object[,se,drop=FALSE],2,
                                  function(x) c(mean(x,na.rm=TRUE))))
        est <- rbind(est,"SE/SD"=est["SE",]/est["SD",])

    }
    if (!is.null(confint)) {
        if (length(confint)!=2*length(estimate)) stop("'confint' should be of length 2*length(estimate).")
        Coverage <- c()
        for (i in seq_along(estimate)) {
            Coverage <- c(Coverage,
                          mean((object[,confint[2*(i-1)+1]]<true[i]) & (object[,confint[2*i]]>true[i]),na.rm=TRUE))
        }
        est <- rbind(est,Coverage=Coverage)
    }
    if (!is.null(names)) colnames(est) <- names
    return(est)
}
