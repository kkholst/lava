uniform.lvm <- function(n,mu,var,...) (mu+(runif(n,-1,1)*sqrt(12)/2*sqrt(var)))
uniform1.lvm <- function(n,...) runif(n,0,1)
binomial.lvm <- function(n,mu,...) rbinom(n,1,tigol(mu))
logit.lvm <- binomial.lvm
probit.lvm <- function(n,mu,var=1,...) rbinom(n,1,pnorm(mu,sd=sqrt(var)))
##probit.lvm <- function(n,mu=0,var=1,...) (rnorm(n,mu,sqrt(var))>0)*1
normal.lvm <- function(n,mu=0,var=1) rnorm(n,mu,sqrt(var))
poisson.lvm <- function(n,mu,...) rpois(n,exp(mu))

"sim" <- function(x,...) UseMethod("sim")
 
sim.lvmfit <- function(x,n=nrow(model.frame(x)),p=pars(x),normal=TRUE,cond=TRUE,...) {
  m <- Model(x)
  if (cond & (nrow(model.frame(x))==n)) {
    X <- exogenous(x)
    mydata <- model.frame(x)
    for (pred in X) {
      distribution(m, pred) <- list(mydata[,pred])
    }
  }
  sim(m,n=n,p=p,normal=normal,cond=cond)
}

sim.lvm <- function(x,n=100,p=NULL,normal=FALSE,cond=FALSE,sigma=1,rho=.5,...) {
  require("MASS")
  if (length(constrain(x))>0) warning("Not using 'constrain' specifications in simulation'")
  x$constrain <- NULL
  
  nn <- vars(x)
  k <- length(nn)
  mu <- unlist(lapply(x$mean, function(l) ifelse(is.na(l)|is.character(l),0,l)))
  xf <- intersect(unique(parlabels(x)),exogenous(x))
  xfix <- c(randomslope(x),xf); if (length(xfix)>0) normal <- FALSE
    
  if (!is.null(p)) {
    if (length(p)!=(index(x)$npar+index(x)$npar.mean)) {
      p0 <- p
      p <- rep(1, index(x)$npar+index(x)$npar.mean)
      p[1:index(x)$npar.mean] <- sigma
      p[offdiags(x)] <- rho
      idx1 <- match(names(p0),coef(x,mean=TRUE))
      idx2 <- which(names(p0)%in%coef(x,mean=TRUE))
      p[idx1] <- p0[idx2]
    }
    M <- modelVar(x,p)
    A <- M$A; P <- M$P ##Sigma <- M$P
    if (!is.null(M$v)) mu <- M$v
  } else {
    P <- index(x)$P
    Pfree <- index(x)$P0
    diag(P)[diag(Pfree)==1] <- sigma
    diag(Pfree) <- 0
    P[Pfree==1] <- rho    
    A <- index(x)$M0
    for (i in 1:k)
      for (j in 1:k)
        if (!is.na(x$fix[i,j])) A[i,j] <- x$fix[i,j]      
    ## Sigma <- diag(k)*sigma
    ## for (i in 1:k)
    ##   for (j in i:k) {
    ##     if (x$cov[i,j]==1 & i!=j)          
    ##       Sigma[i,j] <- Sigma[j,i] <- rho
    ##     if (!is.na(x$covfix[i,j]))
    ##       Sigma[i,j] <- Sigma[j,i] <- x$covfix[i,j]
    ##   }    
  }
 
  E <- mvrnorm(n,rep(0,ncol(P)),P) ## Error term for conditional normal distributed variables
  
  ## Simulate exogenous variables (covariates)
  res <- matrix(0,ncol=length(nn),nrow=n); colnames(res) <- nn
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
      yy <- Ey.x + mvrnorm(n,mu=rep(0,ncol(Vy.x)),Sigma=Vy.x)
      res <- cbind(yy, res[,X]); colnames(res) <- c(colnames(Vy.x),X)
      return(res)
    }
    ## Simulate from sim. distribution (Y,X) (mv-normal)
    I <- diag(k)
    IAi <- solve(I-t(A))
    dd <- mvrnorm(n,mu,P)
    res <- dd%*%t(IAi)
    return(data.frame(res))       
  }

  while (length(simuled)<length(vars(x))) {
    leftovers <- setdiff(vars(x),simuled)
    for (i in leftovers) {
      pos <- match(i,vars(x))
      relations <- colnames(A)[A[,pos]!=0]
      if (all(relations%in%simuled)) { ## Only depending on already simulated variables
        ##        mu.i <- 0
        mu.i <- mu[pos]
        for (From in relations) {
          f <- functional(x,i,From)[[1]]
          if (!is.function(f))
            f <- function(x) x
          reglab <- regfix(x)$labels[From,pos]
          if (reglab%in%xfix) {
            mu.i <- mu.i + res[,reglab]*f(res[,From])
          } else {
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
  
  return(data.frame(res))
}
