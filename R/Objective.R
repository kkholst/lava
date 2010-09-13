###{{{ gaussian

`gaussian_objective.lvm` <-
  function(x,p,data,S,mu,n,...) {
    mp <- modelVar(x,p=p,data=data,...)
    C <- mp$C ## Model specific covariance matrix
    xi <- mp$xi ## Model specific mean-vector 
    iC <- Inverse(C,det=TRUE)
    detC <- attributes(iC)$det
    if (n<2) {
      z <- as.numeric(data-xi)
      val <- log(detC) + (t(z)%*%iC%*%z)[1]
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


`gaussian_hessian.lvm` <- function(p,model,n,opt,mu,...) {
  dots <- list(...); dots$weight <- NULL  
  do.call("information", c(list(x=model,p=p,n=n),dots))
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
  ## If not already done, calculate some relevant zero-one matrices and matrix derivatives

  if (!is.null(data)) {  
    if ((nrow(data)<2 | !is.null(weight))| indiv)
    {
      ##pp <- modelPar(x,p)    
      mp <- modelVar(x,p,data=data[1,])
      iC <- Inverse(mp$C,0,det=FALSE)
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
      for (i in 1:NROW(data)) {
        z <- as.numeric(data[i,])
        u <- z-mp$xi
        if (!is.null(weight)) {
          W <- diag(as.numeric(weight[i,]))
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
  iC <- Inverse(C,0,det=FALSE)
  Debug("Sufficient stats.",debug)
  if (!is.null(mu)) {
    W <- tcrossprod(mu-xi)
    T <- S+W
  } else {
    T <- S
  }
  D <- deriv(x, meanpar=attributes(mp)$meanpar, mom=mp, p=p, mu=mu, mean=mean) ##, all=length(constrain(x))>0)
  vec.iC <- as.vector(iC)  
  Grad <- n/2*crossprod(D$dS, as.vector(iC%*%T%*%iC)-vec.iC)
  if (!is.null(mu)) # & mp$npar.mean>0)
    Grad <- Grad - n/2*crossprod(D$dT,vec.iC)
  return(rbind(as.numeric(Grad)))
}


###}}} gaussian

###{{{ gaussian variants

## Maximum Likelihood with numerical gradient + hessian

gaussian0_objective.lvm <- gaussian_objective.lvm

gaussian2_method.lvm <- "NR"
gaussian2_objective.lvm <- gaussian_objective.lvm
gaussian2_gradient.lvm <- gaussian_gradient.lvm
gaussian2_hessian.lvm <- function(p,model,n,data,...) {
  S <- -score(model,p=p,n=n,data=data,...)
  I <- t(S)%*%S
  attributes(I)$grad <- colSums(S)
  return(I)  
}

gaussian3_objective.lvm <- gaussian_objective.lvm
gaussian3_gradient.lvm <- gaussian_gradient.lvm
gaussian3_hessian.lvm <- function(p,model,n,data,...) {
  I <- information(x=model,p=p,n=n,...)
  S <- score(model,p=p,n=n,data=data)
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
`weighted_hessian.lvm` <- function(p,model,n,opt,weight=NULL,data,...) {
  S <- score(model,p=p,n=n,data=data,weight=weight)
  I <- t(S)%*%S
  attributes(I)$grad <- -colSums(S)
  return(I)  
}

weighted_gradient.lvm <-  function(x,p,data,
                                   S,mu,n,weight=NULL,...) {
  val <- -score(x,p=p,S=S,mu=mu,n=n,data=data,weight=weight)
  if (!is.null(nrow(val))) {
    val <- colSums(val)
  }
  val
}

`weighted_objective.lvm` <- function(x,p,data,
                                   S,mu,n,weight=NULL,...) {
  S <- -score(x,p=p,S=S,mu=mu,n=n,data=data,weight=weight)
  S <- matrix(S,ncol=1)
  res <- t(S)%*%S
  return(res[1])
}


weighted2_method.lvm <- "NR"
`weighted2_hessian.lvm` <- function(p,model,n,opt,weight=NULL,data,...) {
  require(numDeriv)
  myfun <- function(p0) weighted2_gradient.lvm(model,p=p0,n=n,data=data,weight=weight)
  jacobian(myfun,p)
}

weighted2_gradient.lvm <-  weighted_gradient.lvm


`weighted2_objective.lvm` <- function(x,p,data,
                                   S,mu,n,weight=NULL,...) {
  S <- -score(x,p=p,S=S,mu=mu,n=n,data=data,weight=weight)
  S <- matrix(S,ncol=1)
  res <- t(S)%*%S
  return(res[1])
}

###}}} Weighted

###{{{ IPW/MSM

###}}}

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
    A <- m.$A
    P <- m.$P
    J <- m.$J
    IAi <- m.$IAi
    npar.reg <- m.$npar.reg; npar <- m.$npar
    G <- J%*%IAi
    detC <- det(C)
    iC <- try(solve(C), silent=TRUE)
    if (detC<0 | inherits(iC, "try-error"))
      return(.Machine$double.xmax)    
    res <- n/2*(log(detC) + tr(S%*%iC) - log(det(S)) - npar)
    res
  }
###}}} ObjectiveSimple

