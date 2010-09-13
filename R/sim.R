normal.lvm <- function(mean,sd,log=FALSE,...) {
  rnormal <- if(log) rlnorm else rnorm
  if (!missing(mean) & !missing(sd)) 
    f <- function(n,...) rnormal(n,mean,sd)
  else
    f <- function(n,mu=0,var=1,...) rnormal(n,mu,sqrt(var))
  return(f)
}
poisson.lvm <- function(lambda) {
 if (!missing(lambda))
    f <- function(n,...) rpois(n,lambda)
 else
   f <- function(n,mu,...) rpois(n,exp(mu))
 return(f)
}
binomial.lvm <- function(link="logit",p) {
  if (!missing(p))
    f <- function(n,...) rbinom(n,1,p)
  else {
    f <- switch(link,
                logit = 
                function(n,mu,...) rbinom(n,1,tigol(mu)),
                cloglog =
                function(n,mu,...) rbinom(n,1,1-exp(-exp(1-mu))),
                function(n,mu,var=1,...) rbinom(n,1,pnorm(mu,sd=sqrt(var)))
                ### function(n,mu=0,var=1,...) (rnorm(n,mu,sqrt(var))>0)*1
                )
  }
  return(f)
}
uniform.lvm <- function(a,b) {
  if (!missing(a) & !missing(b)) 
    f <- function(n,...) runif(n,a,b)
  else
    f <- function(n,mu,var,...)
      (mu+(runif(n,-1,1)*sqrt(12)/2*sqrt(var)))
  return(f)
}
weibull.lvm <- function(scale=1.25,shape=2,cens=Inf,breakties=0) {
  require(survival)
  lambda <- 1/scale
  f <- function(n,mu,var,...) {
    a0 <- function(t) lambda*shape*(lambda*t)^(shape-1)
    A0 <- function(t) (lambda*t)^shape
    A0i <- function(eta) eta^(1/shape)/lambda
    U <- rexp(n, 1) #give everyone a random death time, on the CH scale
    Z <- U*exp(-mu)
    T <- A0i(Z)
    if (breakties!=0)
      T <- T+runif(n,0,breakties)
    if (is.function(cens))
      cens <- cens(n,...)
    Delta <- (T<cens)
    T[!Delta] <- cens[!Delta]
    S <- Surv(T,Delta*1)
    return(S)
  }
  return(f)
}
logit.lvm <- binomial.lvm("logit")
probit.lvm <- binomial.lvm("probit")


"sim" <- function(x,...) UseMethod("sim")
 
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


sim.lvm <- function(x,n=100,p=NULL,normal=FALSE,cond=FALSE,sigma=1,rho=.5,...) {
  require("mvtnorm")
  index(x) <- reindex(x)
  nn <- setdiff(vars(x),parameter(x))
  mu <- unlist(lapply(x$mean, function(l) ifelse(is.na(l)|is.character(l),0,l)))
  xf <- intersect(unique(parlabels(x)),exogenous(x))
  xfix <- c(randomslope(x),xf); if (length(xfix)>0) normal <- FALSE
  if (length(p)!=(index(x)$npar+index(x)$npar.mean)) {
    p0 <- p
    p <- rep(1, index(x)$npar+index(x)$npar.mean)
    p[1:index(x)$npar.mean] <- 0
    p[index(x)$npar.mean + variances(x)] <- sigma
    p[index(x)$npar.mean + offdiags(x)] <- rho
    idx1 <- match(names(p0),coef(x,mean=TRUE,fix=FALSE))
    idx11 <- match(names(p0),coef(x,mean=TRUE,fix=FALSE,labels=TRUE))
    idx2 <- which(names(p0)%in%coef(x,mean=TRUE,fix=FALSE))
    idx22 <- which(names(p0)%in%coef(x,mean=TRUE,fix=FALSE,labels=TRUE))
##    browser()
    if (length(idx1)>0 && !is.na(idx1))      
      p[idx1] <- p0[idx2]
    if (length(idx11)>0 && !is.na(idx11))
      p[idx11] <- p0[idx22]
  }
  M <- modelVar(x,p,data=NULL)
  A <- M$A; P <- M$P ##Sigma <- M$P
  if (!is.null(M$v)) mu <- M$v
 
  E <- rmvnorm(n,rep(0,ncol(P)),P) ## Error term for conditional normal distributed variables
     
  ## Simulate exogenous variables (covariates)
  res <- matrix(0,ncol=length(nn),nrow=n); colnames(res) <- nn
  res <- as.data.frame(res)
  X <- unique(c(exogenous(x, latent=TRUE, index=FALSE),xfix))
  X.idx <- match(X,vars(x))
  if (!is.null(X) && length(X)>0)
  for (i in 1:length(X)) {
    mu.x <- mu[X.idx[i]]
    dist.x <- distribution(x,X[i])[[1]]
    if (is.function(dist.x)) {
      res[,X.idx[i]] <- dist.x(n=n,mu=mu.x,var=P[X.idx[i],X.idx[i]])
    } else {
      if (is.null(dist.x) || is.na(dist.x)) {
        ##        res[,X.idx[i]] <- rnorm(n,mu.x,sd=Sigma[X.idx[i],X.idx[i]]^0.5)
        res[,X.idx[i]] <- mu.x+E[,X.idx[i]]
      } else {
        res[,X.idx[i]] <- dist.x ## Deterministic
      }
    }
  }
  simuled <- X

  if ( normal | ( is.null(distribution(x)) & is.null(functional(x)) ) ) { ## || all(is.na(distribution(x))) ) {
    if(cond) { ## Simulate from conditional distribution of Y given X
      mypar <- pars(x,A,P,mu)
      pp <- predict(x, mypar, data.frame(res))
      Ey.x <- t(attributes(pp)$Ey.x)
      Vy.x <- attributes(pp)$cond.var
      yy <- Ey.x + rmvnorm(n,mean=rep(0,ncol(Vy.x)),sigma=Vy.x)
      res <- cbind(yy, res[,X]); colnames(res) <- c(colnames(Vy.x),X)
      return(res)
    }
    ## Simulate from sim. distribution (Y,X) (mv-normal)
    I <- diag(length(nn))
    IAi <- solve(I-t(A))
    dd <- rmvnorm(n,mu,P)
    res <- dd%*%t(IAi)
    return(data.frame(res))       
  }

  
  
  xconstrain.idx <- unlist(lapply(lapply(constrain(x),function(z) attributes(z)$args),function(z) length(intersect(z,index(x)$manifest))>0))  
  xconstrain <- intersect(unlist(lapply(constrain(x),function(z) attributes(z)$args)),index(x)$manifest)
  if (!all(xconstrain %in% index(x)$exogenous)) stop("Non-linear constraint only allowed via covariates")
  if (length(xconstrain>0))
  for (i in which(xconstrain.idx)) {
    ff <- constrain(x)[[i]]
    myargs <- attributes(ff)$args
    D <- matrix(0,n,length(myargs))
    for (j in 1:ncol(D)) {
      if (myargs[j]%in%xconstrain)
        D[,j] <- res[,myargs[j]]
      else
        D[,j] <- M$parval[[myargs[j]]]
    }
    res[,names(xconstrain.idx)[i]] <- apply(D,1,ff)
  }
  xconstrain.par <- names(xconstrain.idx)[xconstrain.idx]  
  covparnames <- unique(as.vector(covariance(x)$labels))  
  if (any(xconstrain.par%in%covparnames)) {
    mu0 <- rep(0,ncol(P))
    P0 <- P
    E <- t(sapply(1:n,function(idx) {
      for (i in intersect(xconstrain.par,covparnames)) {
        P0[covariance(x)$labels==i] <- res[idx,i]
      }
      return(rmvnorm(1,mu0,P0))
    }))
  }


    while (length(simuled)<length(nn)) {
    leftovers <- setdiff(nn,simuled)
    
    for (i in leftovers) {
      pos <- match(i,vars(x))
      relations <- colnames(A)[A[,pos]!=0]
      if (all(relations%in%simuled)) { ## Only depending on already simulated variables
        ##        mu.i <- 0
        if (x$mean[[pos]]%in%xconstrain.par) {
          mu.i <- res[,x$mean[[pos]] ]
        } else {
          mu.i <- mu[pos]
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
        if (!is.function(dist.i))
          res[,pos] <- mu.i + E[,pos]
##          res[,pos] <- rnorm(n,mu.i,sd=Sigma[pos,pos]^0.5)
        else {
          res[,pos] <- dist.i(n=n,mu=mu.i,var=P[pos,pos])
        }
        simuled <- c(simuled,i)
      }
    }
  }
  res <- res[,nn,drop=FALSE]

  myhooks <- gethook("sim.hooks")
  count <- 3
  for (f in myhooks) {
    count <- count+1
    res <- do.call(f, list(x=x,data=res))
  }         
  
  return(data.frame(res))
}
