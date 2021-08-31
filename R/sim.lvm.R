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
##' multinomial.lvm
##' beta.lvm
##' normal.lvm mvn.lvm
##' lognormal.lvm
##' gaussian.lvm
##' GM2.lvm
##' GM3.lvm
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
##' ones.lvm Binary.lvm binary.lvm
##' Sequence.lvm
##' none.lvm constant.lvm id.lvm
##' @usage
##' \method{sim}{lvm}(x, n = NULL, p = NULL, normal = FALSE, cond = FALSE,
##' sigma = 1, rho = 0.5, X = NULL, unlink=FALSE, latent=TRUE,
##' use.labels = TRUE, seed=NULL, ...)
##' @param x Model object
##' @param n Number of simulated values/individuals
##' @param p Parameter value (optional)
##' @param normal Logical indicating whether to simulate data from a
##' multivariate normal distribution conditional on exogenous variables hence
##' ignoring functional/distribution definition
##' @param cond for internal use
##' @param sigma Default residual variance (1)
##' @param rho Default covariance parameter (0.5)
##' @param X Optional matrix of fixed values of variables (manipulation)
##' @param unlink Return Inverse link transformed data
##' @param latent Include latent variables (default TRUE)
##' @param use.labels convert categorical variables to factors before applying transformation
##' @param seed Random seed
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
##' ### Beta
##' ##################################################
##' m <- lvm()
##' distribution(m,~y) <- beta.lvm(alpha=2,beta=1)
##' var(sim(m,100,"y,y"=2))
##' distribution(m,~y) <- beta.lvm(alpha=2,beta=1,scale=FALSE)
##' var(sim(m,100))
##'
##' ##################################################
##' ### Transform
##' ##################################################
##' m <- lvm()
##' transform(m,xz~x+z) <- function(x) x[1]*(x[2]>0)
##' regression(m) <- y~x+z+xz
##' d <- sim(m,1e3)
##' summary(lm(y~x+z + x*I(z>0),d))
##'
##' ##################################################
##' ### Non-random variables
##' ##################################################
##' m <- lvm()
##' distribution(m,~x+z+v+w) <- list(Sequence.lvm(0,5),## Seq. 0 to 5 by 1/n
##'                                Binary.lvm(),       ## Vector of ones
##'                                Binary.lvm(0.5),    ##  0.5n 0, 0.5n 1
##'                                Binary.lvm(interval=list(c(0.3,0.5),c(0.8,1))))
##' sim(m,10)
##'
##' ##################################################
##' ### Cox model
##' ### piecewise constant hazard
##' ################################################
##' m <- lvm(t~x)
##' rates <- c(1,0.5); cuts <- c(0,5)
##' ## Constant rate: 1 in [0,5), 0.5 in [5,Inf)
##' distribution(m,~t) <- coxExponential.lvm(rate=rates,timecut=cuts)
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
##' ##################################################
##' ### Cox model
##' ### piecewise constant hazard, gamma frailty
##' ##################################################
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
##' ## Equivalent via transform (here with Aalens additive hazard model)
##' m <- lvm(y~x)
##' distribution(m,~y) <- aalenExponential.lvm(rate=rates,timecut=cuts)
##' distribution(m,~z) <- Gamma.lvm(rate=1,shape=1)
##' transform(m,t~y+z) <- prod
##' sim(m,10)
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
##' distribution(m,~y1+y2+y3,TRUE) <- function(n,...) rmvn0(n,sigma=diag(3)+1)
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
##' ### Regression design in other parameters
##' ##################################################
##' ## Variance heterogeneity
##' m <- lvm(y~x)
##' distribution(m,~y) <- function(n,mean,x) rnorm(n,mean,exp(x)^.5)
##' if (interactive()) plot(y~x,sim(m,1e3))
##' ## Alternaively, calculate the standard error directly
##' addvar(m) <- ~sd ## If 'sd' should be part of the resulting data.frame
##' constrain(m,sd~x) <- function(x) exp(x)^.5
##' distribution(m,~y) <- function(n,mean,sd) rnorm(n,mean,sd)
##' if (interactive()) plot(y~x,sim(m,1e3))
##'
##' ## Regression on variance parameter
##' m <- lvm()
##' regression(m) <- y~x
##' regression(m) <- v~x
##' ##distribution(m,~v) <- 0 # No stochastic term
##' ## Alternative:
##' ## regression(m) <- v[NA:0]~x
##' distribution(m,~y) <- function(n,mean,v) rnorm(n,mean,exp(v)^.5)
##' if (interactive()) plot(y~x,sim(m,1e3))
##'
##' ## Regression on shape parameter in Weibull model
##' m <- lvm()
##' regression(m) <- y ~ z+v
##' regression(m) <- s ~ exp(0.6*x-0.5*z)
##' distribution(m,~x+z) <- binomial.lvm()
##' distribution(m,~cens) <- coxWeibull.lvm(scale=1)
##' distribution(m,~y) <- coxWeibull.lvm(scale=0.1,shape=~s)
##' eventTime(m) <- time ~ min(y=1,cens=0)
##'
##' if (interactive()) {
##'     d <- sim(m,1e3)
##'     require(survival)
##'     (cc <- coxph(Surv(time,status)~v+strata(x,z),data=d))
##'     plot(survfit(cc) ,col=1:4,mark.time=FALSE)
##' }
##'
##' ##################################################
##' ### Categorical predictor
##' ##################################################
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
##'
##' transform(m,v2~v) <- function(x) x=='A'
##' sim(m,10)
##'
##' ##################################################
##' ### Pre-calculate object
##' ##################################################
##' m <- lvm(y~x)
##' m2 <- sim(m,'y~x'=2)
##' sim(m,10,'y~x'=2)
##' sim(m2,10) ## Faster
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
sim.lvm <- function(x,n=NULL,p=NULL,normal=FALSE,cond=FALSE,sigma=1,rho=.5,
            X=NULL,unlink=FALSE,latent=TRUE,use.labels=TRUE,seed=NULL,...) {

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


    v.env <- c("A","M","P","PP","PPdiag","xx","vv","mdist","mdistnam","mii",
              "nn","mu","xf","xfix","X",
              "vartrans","multitrans","multitrans.idx",
              "X.idx","ii.mvn","xconstrain.idx","xconstrain",
              "xconstrain.par","covparnames","exo_constrainY")

    setup <- is.null(n) && is.null(X) ## Save environment (variables v.env) and return sim object
    loadconfig <- !is.null(x$sim.env) && !setup && (length(list(...))==0 && length(p)==0)
    Yfix <- NULL

    if (loadconfig) {
        for (v in setdiff(v.env,"X")) assign(v, x$sim.env[[v]])
        if (is.null(X)) X <-  x$sim.env[['X']]
    } else {
        if (!is.null(n) && n<1) return(NULL)
        p <- c(p,unlist(list(...)))
        xx <- exogenous(x)

        if (!is.null(X)) {
            if (is.null(n))
                n <- nrow(X)
            if (!is.null(colnames(X))) {
                yfix <- setdiff(colnames(X),xx)
                if (length(yfix)>0) Yfix <- X[,yfix,drop=FALSE]
                xx0 <- intersect(xx,colnames(X))
                if (length(xx0)>0)
                    X <- as.matrix(X[,xx0,drop=FALSE])
            }
        } else {
            if (!is.null(p)) {
                i1 <- unique(na.omit(c(match(names(p),xx),
                                       match(names(p),paste0(xx,lava.options()$symbol[2],xx)))))
                covariance(x) <- xx[i1]
            }
        }
        ##  index(x) <- reindex(x)
        vv <- vars(x)
        nn <- setdiff(vv,parameter(x))
        mu <- unlist(lapply(x$mean, function(l) ifelse(is.na(l)|is.character(l),0,l)))
        xf <- intersect(unique(parlabels(x)),xx)
        xfix <- c(randomslope(x),xf); if (length(xfix)>0) normal <- FALSE
        ## Match parameter names
        if ((!is.null(names(p)) && all(!is.na(names(p)))) || length(p)!=(index(x)$npar+index(x)$npar.mean+index(x)$npar.ex) | is.null(names(p))) {
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
                if (length(idx1)>0 && !any(is.na(idx1)))
                    p[idx1] <- p0[idx2]
                if (length(idx11)>0 && !any(is.na(idx11)))
                    p[idx11] <- p0[idx22]
            }
        }
        M <- modelVar(x,p,data=NULL)
        A <- M$A; P <- M$P
        if (!is.null(M$v)) mu <- M$v

        ## Square root of residual variance matrix
        PP <- with(svd(P), v%*%diag(sqrt(d),nrow=length(d))%*%t(u))
        ## Multivariate distributions
        mdist <- distribution(x,multivariate=TRUE)$var
        mdistnam <- names(mdist)
        mii <- match(mdistnam,vars(x))
        if (length(distribution(x))>0 ) {
            ii <- match(names(distribution(x)),vv)
            ii.mvn <- setdiff(seq(ncol(P)),c(ii,mii))
        } else {
            ii.mvn <- seq(ncol(P))
        }
        PPdiag <- sum(abs(offdiag(PP[ii.mvn,ii.mvn,drop=FALSE])^2))<1e-20
    }

    if (!setup) {
        E <- matrix(0,ncol=ncol(P),nrow=n)
        if (length(ii.mvn)>0) {
            ## Error term for conditional normal distributed variables
            if (PPdiag) {
                for (i in ii.mvn) E[,i] <- rnorm(n,sd=PP[i,i])
            } else {
                E[,ii.mvn] <-  matrix(rnorm(length(ii.mvn)*n),ncol=length(ii.mvn))%*%PP[ii.mvn,ii.mvn,drop=FALSE]
            }
        }

        if (length(mdistnam)>0) {
            fun <- distribution(x,multivariate=TRUE)$fun
            for (i in seq_along(fun)) {
                mv <- names(which(unlist(mdist)==i))
                ii <- match(mv,vv)
                E[,ii] <- distribution(x,multivariate=TRUE)$fun[[i]](n,p=p,object=x) # ,...)
            }
        }

        ## Simulate exogenous variables (covariates)
        res <- matrix(0,ncol=length(nn),nrow=n)
        colnames(res) <- nn
    }

    if (!loadconfig) {
        vartrans <- names(x$attributes$transform)
        multitrans <- multitrans.idx <- NULL
        if (length(x$attributes$multitransform)>0) {
            multitrans <- unlist(lapply(x$attributes$multitransform,function(z) z$y))
            for (i in (seq_along(x$attributes$multitransform))) {
                multitrans.idx <- c(multitrans.idx,rep(i,length(x$attributes$multitransform[[i]]$y)))
            }
        }
        xx <- unique(c(exogenous(x, latent=FALSE, index=TRUE),xfix))
        xx <- setdiff(xx,vartrans)

        X.idx <- match(xx,vv)
    }

    if (!setup) {
        res[,X.idx] <- t(mu[X.idx]+t(E[,X.idx]))
        if (is.null(X) || NCOL(X)<length(xx)) {
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
        }
        if (!is.null(X)) {
            ii <- match(colnames(X),vv)
            res0 <- res
            for (i in seq(ncol(X)))
                res[,ii[i]] <- X[,i]
        }
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
        I <- diag(nrow=length(nn))
        IAi <- Inverse(I-t(A))
        colnames(E) <- vv
        dd <- t(apply(heavytail.sim.hook(x,E),1,function(x) x+mu))
        res <- dd%*%t(IAi)
        colnames(res) <- vv
    } else {

        if (!loadconfig) {
            xc <- index(x)$vars
            xconstrain.idx <- unlist(lapply(lapply(constrain(x),function(z) attributes(z)$args),function(z) length(intersect(z,xc))>0))
            xconstrain <- intersect(unlist(lapply(constrain(x),function(z) attributes(z)$args)),xc)
            xconstrain.par <- names(xconstrain.idx)[xconstrain.idx]
            covparnames <- unique(as.vector(covariance(x)$labels))
            exo_constrainY <- intersect(exogenous(x),names(x$constrainY))
        }

        if (setup) {
            sim.env <- c()
            sim.env[v.env] <- list(NULL)
            for (v in v.env) if (!is.null(get(v))) sim.env[[v]] <- get(v)
            x$sim.env <- sim.env
            return(x)
        }

      if (length(xconstrain)>0) {
        for (i in which(xconstrain.idx)) {
          if (names(xconstrain.idx[i]) %in% nn) { ## A parameter and not a variable
            ff <- constrain(x)[[i]]
            myargs <- attributes(ff)$args
            D <- matrix(0,n,length(myargs))
            for (j in seq_len(ncol(D))) {
              if (myargs[j]%in%xconstrain)
                D[,j] <- res[,myargs[j]]
              else
                D[,j] <- M$parval[[myargs[j]]]
            }
            val <- try(val <- ff(D), silent=TRUE)
            if (inherits(val,"try-error") || NROW(val)<n) apply(D,1,ff)
            res <- cbind(res, val); colnames(res)[ncol(res)] <- names(xconstrain.idx)[i]
          }
        }
      }

        if (any(xconstrain.par%in%covparnames)) {
            mu0 <- rep(0,ncol(P))
            P0 <- P
            E <- t(sapply(seq_len(n),function(idx) {
              for (i in intersect(xconstrain.par,covparnames)) {
                  P0[covariance(x)$labels==i] <- res[idx,i]
                }
                PP <- with(svd(P0), v%*%diag(sqrt(d),nrow=length(d))%*%t(u))
                return(mu0+rbind(rnorm(ncol(P0)))%*%PP)
            }))
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
        yconstrain <- unlist(lapply(xconstrain,function(x) x$endo))
      for (i in exo_constrainY) {
            cc <- x$constrainY[[i]]
            args <- cc$args
            args <- if (is.null(args) || length(args)==0) res[,i] else res[,args]
            res[,i] <- cc$fun(args,p) # ,...)
        }

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
        itercount <- 0
        while (length(simuled)<length(nn)) {
            leftoversPrev <- leftovers
            leftovers <- setdiff(nn,simuled)
            if (!is.null(leftoversPrev) && length(leftoversPrev)==length(leftovers)) {
                if (itercount>0)
                stop("Infinite loop (feedback).")
                itercount <- itercount+1
            }
            for (i in leftovers) {
                if (length(Yfix)>0 && i %in% colnames(Yfix)) {
                    if (NROW(Yfix) == 1) {
                        res[,i] <- rep(Yfix[,i],length.out=nrow(res))
                    } else                     
                        res[,i] <- rep(Yfix[,i,drop=TRUE],length.out=nrow(res))
                    simuled <- c(simuled,i)
                    next
                }
                if (i%in%vartrans) {
                    xtrans <- x$attributes$transform[[i]]$x
                    if (all(xtrans%in%c(simuled,names(parvals))))  {
                        xtr <- res[,xtrans,drop=FALSE]
                        if (use.labels) {
                            lb <- x$attributes$labels
                            lb.idx <- na.omit(match(names(lb),xtrans))
                            ## For categorical variables turn them into factors so we can
                            ## use the actual labels in function calls/transform
                            if (length(lb.idx)>0) {
                                xtr <- as.data.frame(xtr)
                                for (lb0 in lb.idx) {
                                    lab <- lb[[names(xtr)[lb0]]]
                                    xtr[,lb0] <- factor(xtr[,lb0],levels=seq_along(lab)-1,labels=lab)
                                }
                            }
                        }
                        suppressWarnings(yy <- with(x$attributes$transform[[i]], fun(xtr))) ##fun(res[,xtrans])))
                        if (NROW(yy) != NROW(res)) { ## apply row-wise
                            res[,i] <- with(x$attributes$transform[[i]], ##apply(res[,xtrans,drop=FALSE],1,fun))
                                           apply(xtr,1,fun))
                        } else {
                            colnames(yy) <- NULL
                            res[,i] <- yy
                        }
                        simuled <- c(simuled,i)
                    }
                } else if (i%in%multitrans) {
                    idx0 <- match(i,multitrans)
                    idx <- multitrans.idx[idx0]
                    mtval <- x$attributes$multitransform[[idx]]
                    if (all(mtval$x%in%simuled)) {
                        res[,mtval$y] <- mtval$fun(res[,mtval$x])
                        simuled <- c(simuled,mtval$y)
                        break;
                    }
                } else {
                    ipos <- which(yconstrain%in%i)
                    if (length(ipos)==0 || all(xconstrain[[ipos]]$exo%in%simuled)) {
                        pos <- match(i,vv)
                        relations <- colnames(A)[A[,pos]!=0]
                        simvars <- x$attributes$simvar[[i]]
                        dist.i <- distribution(x,i)[[1]] ## User-specified distribution function
                        dist.xx <- NULL
                        if (is.function(dist.i)) {
                            dist.args0 <- names(formals(dist.i))
                            dist.args <- setdiff(dist.args0,c("n","mean","mu","var","..."))
                            dist.xx <- intersect(names(res),dist.args) ## Variables influencing distribution
                        }
                        if (all(c(relations,simvars,dist.xx)%in%simuled)) { ## Only depending on already simulated variables
                            if (x$mean[[pos]]%in%xconstrain.par && length(ipos)==0) {
                                mu.i <- res[,x$mean[[pos]] ]
                            } else {
                                mu.i <- mu[pos]
                            }
                            if (length(ipos)>0) {
                                pp <- unlist(M$parval[xconstrain[[ipos]]$warg])
                                myidx <- with(xconstrain[[ipos]],order(c(wargidx,exoidx)))
                                ## myidx <- with(xconstrain[[ipos]],
                                ##              match(attr(func,"args"), c(warg,exo)))
                                X <- with(xconstrain[[ipos]],
                                          if (length(pp)>0)
                                              cbind(rbind(pp)%x%cbind(rep(1,nrow(res))),
                                                    res[,exo,drop=FALSE])
                                          else res[,exo,drop=FALSE])
                                yy <- try(with(xconstrain[[ipos]],
                                               func(X[,myidx,drop=FALSE])),silent=TRUE)
                                if (NROW(yy) != NROW(res)) { ## apply row-wise
                                    mu.i <- #mu.i +
                                        with(xconstrain[[ipos]],
                                             apply(res[,exo,drop=FALSE],1,
                                                   function(x) func(
                                                            unlist(c(pp,x))[myidx])))
                                } else {
                                    mu.i <- ##mu.i+
                                        yy
                                }
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
                                new.args <- list(n=n)
                                mu.arg <- intersect(c("mean","mu"),dist.args0)
                                if (length(mu.arg)>0) {
                                    new.args <- c(new.args,list(mu.i))
                                    names(new.args)[length(new.args)] <- mu.arg[1]
                                }
                                var.arg <- intersect(c("var"),dist.args0)
                                if (length(var.arg)>0) {
                                    new.args <- c(new.args,list(P[pos,pos]))
                                    names(new.args)[length(new.args)] <- var.arg[1]
                                }
                                for (jj in dist.xx) {
                                    new.args <- c(new.args,list(res[,jj,drop=TRUE]))
                                    names(new.args)[length(new.args)] <- jj
                                }
                                res[,pos] <- do.call(dist.i,new.args)
                                if (unlink)
                                    resunlink[,pos] <- mu.i
                            }

                            if (length(x$constrainY)>0 && i%in%names(x$constrainY)) {
                                cc <- x$constrainY[[i]]
                                args <- cc$args
                                args <- if (length(args)==0)
                                           res[,pos]
                                       else {
                                           ii <- intersect(names(M$parval),args)
                                           args0 <- args
                                           args <- res[,intersect(args0,colnames(res)),drop=FALSE]

                                           if (length(ii)>0) {
                                               pp <- rbind(unlist(M$parval[ii]))%x%cbind(rep(1,n))
                                               colnames(pp) <- ii
                                               args <- cbind(res,pp)[,args0,drop=FALSE]
                                           }
                                           args
                                       }
                                res[,pos] <- cc$fun(args,p) # ,...)
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
    if (!latent && length(latent(x))>0) return(subset(res[,-which(colnames(res)%in%latent(x))]))
    return(res)
}



##' @export
simulate.lvm <- function(object,nsim,seed=NULL,...) {
    sim(object,nsim,seed=seed,...)
}

##' @export
simulate.lvmfit <- function(object,nsim,seed=NULL,...) {
    sim(object,nsim,seed=seed,...)
}

