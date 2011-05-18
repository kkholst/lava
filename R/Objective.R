###{{{ gaussian

gaussian_method.lvm <- "nlminb2"
`gaussian_objective.lvm` <-
  function(x,p,data,S,mu,n,...) {
    mp <- modelVar(x,p=p,data=data,...)
    C <- mp$C ## Model specific covariance matrix
    xi <- mp$xi ## Model specific mean-vector
    iC <- Inverse(C,tol=1e-7,det=TRUE)
    detC <- attributes(iC)$det
    if (n<2) {
      z <- as.numeric(data-xi)
      val <- log(detC) + tcrossprod(z,crossprod(z,iC))[1]
        ##(t(z)%*%iC%*%z)[1]
      return(0.5*val)      
    }
    if (!is.null(mu)){
      ##      if (is.null(mu)) mu <- xi
      W <- suppressMessages(tcrossprod(mu-xi))
      T <- S+W
    } else {
      T <- S
    }
    res <- n/2*log(detC) + n/2*tr(T%*%iC) ## Objective function (Full Information ML)
    return(res)
  }


`gaussian_hessian.lvm` <- function(x,p,n,mu,opt,...) {
  dots <- list(...); dots$weight <- NULL  
  do.call("information", c(list(x=x,p=p,n=n),dots))
}

gaussian_gradient.lvm <-  function(x,p,data,S,mu,n,...) {
  dots <- list(...); dots$weight <- NULL
  if (n>2) data <- NULL
  val <- -gaussian_score.lvm(x,p=p,S=S,mu=mu,n=n,data=data,reindex=FALSE,...)
  if (!is.null(nrow(val))) {
    val <- colSums(val)
  }
  val
}

gaussian_score.lvm <- function(x, data, p, S, n, mu=NULL, weight=NULL, debug=FALSE, reindex=FALSE, mean=TRUE, constrain=TRUE, indiv=FALSE,...) {

  if (!is.null(data)) {  
    if ((nrow(data)<2 | !is.null(weight))| indiv)
    {
      ##pp <- modelPar(x,p)    
      mp <- modelVar(x,p,data=data[1,])
      iC <- Inverse(mp$C,det=FALSE)
##      D <- with(pp, deriv(x, meanpar=meanpar, p=p, mom=mp, mu=NULL)) ##, all=length(constrain(x))>0))
      MeanPar <- attributes(mp)$meanpar
      D <- with(attributes(mp), deriv(x, meanpar=MeanPar, p=pars, mom=mp, mu=NULL)) ##, all=length(constrain(x))>0))
      ##      D <- with(pp, deriv(x, meanpar=meanpar, mom=mp, mu=NULL))
      Debug("after deriv.", debug)
      myvars <- (index(x)$manifest)
      if (NCOL(data)!=length(myvars)) {
        data <- subset(data,select=myvars)
      }
      score <- matrix(ncol=length(p),nrow=NROW(data))
      score0 <- -1/2*as.vector(iC)%*%D$dS
      if (!is.null(weight)) {
        W0 <- diag(length(myvars))        
        widx <- match(colnames(weight),myvars)
      }
      
      for (i in 1:NROW(data)) {
        z <- as.numeric(data[i,])
        u <- z-mp$xi
        if (!is.null(weight)) {
          W <- W0; diag(W)[widx] <- as.numeric(weight[i,])
          score[i,] <- 
            as.numeric(crossprod(u,iC%*%W)%*%D$dxi +
                       -1/2*(as.vector((iC
                                        - iC %*% tcrossprod(u)
                                        %*% iC)%*%W)) %*% D$dS
                       )
        } else {
          score[i,] <-
            as.numeric(score0 + crossprod(u,iC)%*%D$dxi +
                       1/2*as.vector(iC%*%tcrossprod(u)%*%iC)%*%D$dS)
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
  iC <- Inverse(C,det=FALSE)
  Debug("Sufficient stats.",debug)
  if (!is.null(mu) & !is.null(xi)) {
    W <- tcrossprod(mu-xi)
    T <- S+W
  } else {
    T <- S
  }
  D <- deriv(x, meanpar=attributes(mp)$meanpar, mom=mp, p=p, mu=mu, mean=mean) ##, all=length(constrain(x))>0)
  vec.iC <- as.vector(iC)  
  Grad <- n/2*crossprod(D$dS, as.vector(iC%*%T%*%iC)-vec.iC)
  if (!is.null(mu) & !is.null(xi)) # & mp$npar.mean>0)
    Grad <- Grad - n/2*crossprod(D$dT,vec.iC)
  return(rbind(as.numeric(Grad)))
}

###}}} gaussian

###{{{ gaussian variants

## Maximum Likelihood with numerical gradient + hessian
gaussian0_objective.lvm <- gaussian_objective.lvm

gaussian1_objective.lvm <- gaussian_objective.lvm
gaussian1_gradient.lvm <- function(...) gaussian_gradient.lvm(...)
gaussian1_hessian.lvm <- function(x,p,...) {
  ##  browser()
  myg2 <- function(p1) gaussian_gradient.lvm(x,p=p1,...)
  myg3 <- function(p1) numDeriv::jacobian(myg2,p1)
  ##-numDeriv::jacobian(myg,p)
  myg <- function(p1) gaussian_objective.lvm(x,p=p1,...)
  numDeriv::hessian(myg,p)
}

  NULL ##gaussian_hessian.lvm

## BHHH 
gaussian2_method.lvm <- "NR"
gaussian2_objective.lvm <- gaussian_objective.lvm
gaussian2_gradient.lvm <- gaussian_gradient.lvm
gaussian2_hessian.lvm <- function(x,p,n,data,...) {
  S <- -score(x,p=p,n=n,data=data,...)
  I <- t(S)%*%S
  attributes(I)$grad <- colSums(S)
  return(I)  
}
## Sandwich
gaussian3_objective.lvm <- gaussian_objective.lvm
gaussian3_gradient.lvm <- gaussian_gradient.lvm
gaussian3_hessian.lvm <- function(x,p,n,data,...) {
  I <- information(x=x,p=p,n=n,...)
  S <- score(x=x,p=p,n=n,data=data)
  J <- t(S)%*%S
  return(J%*%solve(I)%*%J)  
}

gaussian4_objective.lvm <- gaussian_objective.lvm
gaussian4_gradient.lvm <- gaussian_gradient.lvm
gaussian4_hessian.lvm <- gaussian_hessian.lvm
gaussian4_variance.lvm <- function(x,p,data) {
  matrix(0,ncol=length(p),nrow=length(p))
}

###}}}

###{{{ Weighted

weighted_method.lvm <- "NR"
weighted_gradient.lvm <- function(x,p,data,weight,indiv=FALSE,...) { 
  myvars <- index(x)$manifest
  if (NCOL(data)!=length(myvars))
    data <- subset(data,select=myvars)  
  score <- matrix(ncol=length(p),nrow=NROW(data))
  myy <- index(x)$endogenous
  myx <- index(x)$exogenous
  mynx <- setdiff(myvars,myx)
  W0 <- diag(length(myy))        
  widx <- match(colnames(weight),myy)
  pp <- modelPar(x,p)
  mp <- moments(x,p=p,conditional=TRUE,data=data[1,])
  iC <- Inverse(mp$C,det=FALSE)
  v <- matrix(0,ncol=length(vars(x)),nrow=NROW(data))
  colnames(v) <- vars(x)
  for (i in mynx) v[,i] <- mp$v[i]
  for (i in myx) v[,i] <- data[,i]
  xi <- t(mp$G%*%t(v))
  u <- as.matrix(data)[,myy]-xi
  D <- deriv(x, meanpar=pp$meanpar, 
             p=pp$p, mom=mp, mu=NULL)
  if (NROW(data)==1) {
    W <- W0; diag(W)[widx] <- as.numeric(weight[i,])
    score[i,] <- 
      as.numeric(crossprod(u,iC%*%W)%*%D$dxi +
                 -1/2*(as.vector((iC
                                  - iC %*% tcrossprod(u)
                                  %*% iC)%*%W)) %*% D$dS)
    return(-score)
  }
  score0 <- -0.5*as.vector(iC)%*%D$dS
  Gdv <- mp$G%*%D$dv
  ##  browser()
  for (i in 1:NROW(data)) {    
    W <- W0; diag(W)[widx] <- as.numeric(weight[i,])
    dxi <-      
      (t(as.numeric(v[i,]))%x%diag(nrow=length(myy)))%*%D$dG + Gdv
    score[i,] <- -0.5*as.vector(iC%*%W)%*%D$dS + 
      as.numeric(crossprod(u[i,],iC%*%W)%*%dxi +
                 1/2*as.vector(iC%*%tcrossprod(u[i,])%*%iC%*%W)%*%D$dS)
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


weighted0_method.lvm <- "estfun"
weighted0_gradient.lvm <- function(...) {
  val <- -gaussian_score.lvm(...)
  colSums(val)
}
weighted0_hessian.lvm <- NULL  

weighted2_method.lvm <- "estfun"
weighted2_gradient.lvm <- function(x,p,data,weight,indiv=FALSE,...) {
  myvars <- index(x)$manifest
  if (NCOL(data)!=length(myvars))
    data <- subset(data,select=myvars)  
  score <- matrix(ncol=length(p),nrow=NROW(data))
  myy <- index(x)$endogenous
  myx <- index(x)$exogenous
  mynx <- setdiff(myvars,myx)
  W0 <- diag(length(myy))        
  widx <- match(colnames(weight),myy)
  pp <- modelPar(x,p)
  for (i in 1:NROW(data)) {
    z <- as.matrix(data[i,myy])
    mp <- moments(x,p=p,conditional=TRUE,data=data[i,])
    ##    browser()
    u <- as.numeric(z-mp$xi[,1])
    iC <- Inverse(mp$C,det=FALSE)
    D <- deriv(x, meanpar=pp$meanpar, 
               p=pp$p, mom=mp, mu=NULL)
    W <- W0; diag(W)[widx] <- as.numeric(weight[i,])
    score[i,] <- -0.5*as.vector(iC%*%W)%*%D$dS + 
      as.numeric(crossprod(u,iC%*%W)%*%D$dxi +
                 1/2*as.vector(iC%*%tcrossprod(u)%*%iC%*%W)%*%D$dS)    
  }
  if (indiv) return(-score)
  colSums(-score)
}
weighted2_hessian.lvm <- NULL



###}}} Weighted



