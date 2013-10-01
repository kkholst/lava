
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
##' student.lvm
##' coxGompertz.lvm
##' coxWeibull.lvm
##' coxExponential.lvm
##' aalenExponential.lvm
##' Gamma.lvm gamma.lvm
##' loggamma.lvm
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
##' 
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
##' ### Cox model
##' ### piecewise constant hazard, gamma frailty
##' ##################################################
##' 
##' m <- lvm(y~x+z)
##' rates <- c(1,0.5); cuts <- c(0,5)
##' distribution(m,~y+z) <- list(coxExponential.lvm(rate=rates,timecut=cuts)
##'                              loggamma.lvm(rate=1,shape=1))
##' \dontrun{
##'     d <- sim(m,2e4,p=c("y~x"=0.1)); d$status <- TRUE
##'     plot(timereg::aalen(Surv(y,status)~x,data=d,resample.iid=0,robust=0),spec=1)
##'     L <- approxfun(c(cuts,max(d$y)),f=1,cumsum(c(0,rates*diff(c(cuts,max(d$y))))),method="linear")
##'     curve(L,0,100,add=TRUE,col="blue")
##' }
##' 
##' ## Equivalent via transform (here with Aalens additive hazard model)
##' m <- lvm(y~x)
##' distribution(m,~y) <- aalenExponential.lvm(rate=rates,timecut=cuts)
##' distribution(m,~z) <- gamma.lvm(rate=1,shape=1)
##' transform(m,t~y+z) <- prod
##' sim(m,10)
##' 
##' 
##' 
##' 
##' 
"sim" <- function(x,...) UseMethod("sim")

##' @S3method sim lvmfit
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

##' @S3method sim lvm
sim.lvm <- function(x,n=100,p=NULL,normal=FALSE,cond=FALSE,sigma=1,rho=.5,
                    X,unlink=FALSE,...) {
  if (!missing(X)) {
    n <- nrow(X)
  }
  xx <- exogenous(x)
  if (!is.null(p)) {
    i1 <- na.omit(c(match(names(p),xx),
                    match(names(p),paste(xx,lava.options()$symbol[2],xx,sep=""))))
    if (length(i1)>0) covariance(x) <- xx[i1]
  }
  index(x) <- reindex(x)
  vv <- vars(x)
  nn <- setdiff(vv,parameter(x))
  mu <- unlist(lapply(x$mean, function(l) ifelse(is.na(l)|is.character(l),0,l)))
  xf <- intersect(unique(parlabels(x)),xx)
  xfix <- c(randomslope(x),xf); if (length(xfix)>0) normal <- FALSE


  if (length(p)!=(index(x)$npar+index(x)$npar.mean) | !is.null(names(p))) {
    nullp <- is.null(p)
    p0 <- p
    p <- rep(1, index(x)$npar+index(x)$npar.mean)
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
  A <- M$A; P <- M$P ##Sigma <- M$P
  if (!is.null(M$v)) mu <- M$v
  
##  E <- rmvnorm(n,rep(0,ncol(P)),P) ## Error term for conditional normal distributed variables
  PP <- with(svd(P), v%*%diag(sqrt(d))%*%t(u))
  E <- matrix(rnorm(ncol(P)*n),ncol=ncol(P))%*%PP
  
  ## Simulate exogenous variables (covariates)
  res <- matrix(0,ncol=length(nn),nrow=n); colnames(res) <- nn
  vartrans <- names(attributes(x)$transform)
  xx <- unique(c(exogenous(x, latent=FALSE, index=TRUE),xfix))
  xx <- setdiff(xx,vartrans)
  
  X.idx <- match(xx,vv)
  res[,X.idx] <- t(mu[X.idx]+t(E[,X.idx]))
  if (missing(X)) {
    if (!is.null(xx) && length(xx)>0)
      for (i in 1:length(xx)) {
        mu.x <- mu[X.idx[i]]
        dist.x <- distribution(x,xx[i])[[1]]
        if (is.list(dist.x)) {
          dist.x <- dist.x[[1]]
          if (length(dist.x)==1) dist.x <- rep(dist.x,n)
        }
        if (is.function(dist.x)) {
          res[,X.idx[i]] <- dist.x(n=n,mu=mu.x,var=P[X.idx[i],X.idx[i]])
        } else {
          if (is.null(dist.x) || is.na(dist.x)) {
          } else {
              if (length(dist.x)!=n) stop("'",vv[X.idx[i]], "' fixed at length ", length(dist.x)," != ",n)
              res[,X.idx[i]] <- dist.x ## Deterministic
          }
        }
      }
  } else {
    res[,X.idx] <- X[,xx]
  }
  simuled <- xx
  resunlink <- NULL
  if (unlink) {
    resunlink <- res
  }
  
  if ( normal | ( is.null(distribution(x)) & is.null(functional(x)) & is.null(constrain(x))) ) { 
    if(cond) { ## Simulate from conditional distribution of Y given X
      mypar <- pars(x,A,P,mu)
      Ey.x <- predict(x, mypar, data.frame(res))
      Vy.x <- attributes(Ey.x)$cond.var
##      yy <- Ey.x + rmvnorm(n,mean=rep(0,ncol(Vy.x)),sigma=Vy.x)
      PP <- with(svd(Vy.x), v%*%diag(sqrt(d))%*%t(u))
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
    ## if (length(xconstrain>0))
    ##   for (i in which(xconstrain.idx)) {
    ##     ff <- constrain(x)[[i]]
    ##     myargs <- attributes(ff)$args
    ##     D <- matrix(0,n,length(myargs))
    ##     for (j in 1:ncol(D)) {
    ##       if (myargs[j]%in%xconstrain)
    ##         D[,j] <- res[,myargs[j]]
    ##       else
    ##         D[,j] <- M$parval[[myargs[j]]]
    ##     }
    ##     res[,names(xconstrain.idx)[i]] <- apply(D,1,ff)
    ##   }
    
    xconstrain.par <- names(xconstrain.idx)[xconstrain.idx]  
    covparnames <- unique(as.vector(covariance(x)$labels))  
    if (any(xconstrain.par%in%covparnames)) {
      mu0 <- rep(0,ncol(P))
      P0 <- P
      E <- t(sapply(1:n,function(idx) {
        for (i in intersect(xconstrain.par,covparnames)) {
          P0[covariance(x)$labels==i] <- res[idx,i]
        }
##        return(rmvnorm(1,mu0,P0))
        PP <- with(svd(P0), v%*%diag(sqrt(d))%*%t(u))
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
    yconstrain <- unlist(lapply(xconstrain,function(x) x$endo))

    res <- data.frame(res)
    if (length(vartrans)>0) {
        parvals <- parpos(x)$parval
        if (length(parvals)>0) {
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

      if (!is.null(leftoversPrev) && length(leftoversPrev)==length(leftovers)) stop("Infinite loop (probably problem with 'transform' call in model: Outcome variable should not affect other variables in the model)")
      for (i in leftovers) {
        if (i%in%vartrans) {
          xtrans <- attributes(x)$transform[[i]]$x
          if (all(xtrans%in%c(simuled,names(parvals))))  {
            suppressWarnings(yy <- with(attributes(x)$transform[[i]],fun(res[,xtrans])))
            if (length(yy) != NROW(res)) { ## apply row-wise
              res[,i] <- with(attributes(x)$transform[[i]],apply(res[,xtrans,drop=FALSE],1,fun))
            } else {
              res[,i] <- yy
            }
            simuled <- c(simuled,i)
          }
        } else {

          ipos <- which(i%in%yconstrain)

          if (length(ipos)==0 || all(xconstrain[[ipos]]$exo%in%simuled)) {
            pos <- match(i,vv)
            relations <- colnames(A)[A[,pos]!=0]
            
            if (all(relations%in%simuled)) { ## Only depending on already simulated variables
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
                f <- function(x) x
              reglab <- regfix(x)$labels[From,pos]
              if (reglab%in%c(xfix,xconstrain.par)) {
                mu.i <- mu.i + res[,reglab]*f(res[,From])
              }
              else {
                mu.i <- mu.i + A[From,pos]*f(res[,From])
              }
            }
            dist.i <- distribution(x,i)[[1]]
            if (!is.function(dist.i)) {
              res[,pos] <- mu.i + E[,pos]
              if (unlink)
                resunlink[,pos] <- res[,pos]
            }
            else {
              res[,pos] <- dist.i(n=n,mu=mu.i,var=P[pos,pos])
              if (unlink)
                resunlink[,pos] <- mu.i
            }
            simuled <- c(simuled,i)
          }
          }
        }
      }
    }
    res <- res[,nn,drop=FALSE]
  }

  ## for (i in seq_len(length(x$constrainY))) {
  ##   cc <- x$constrainY[[i]]
  ##   args <- attributes(x$constrainY[[i]])$args
  ##   nam <- names(x$constrainY)[[i]]
  ##   newcol <- cbind(apply(res[,args,drop=FALSE],1,cc)); colnames(newcol) <- nam
  ##   res <- cbind(res,newcol)
  ## }     

  res <- as.data.frame(res)
  myhooks <- gethook("sim.hooks")
  for (f in myhooks) {
    res <- do.call(f, list(x=x,data=res))
  }         
  if (unlink) res <- resunlink
  return(as.data.frame(res))
}



##' @S3method simulate lvm
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

##' @S3method simulate lvmfit
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

rmvn <- function(n,mu=rep(0,ncol(Sigma)),Sigma=diag(2)+1,...) {
    PP <- with(svd(Sigma), v%*%diag(sqrt(d))%*%t(u))
    res <- matrix(rnorm(ncol(Sigma)*n),ncol=ncol(Sigma))%*%PP
    return(res+cbind(rep(1,n))%*%mu)
}
