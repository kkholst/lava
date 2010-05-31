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
  val <- -score(x,p=p,S=S,mu=mu,n=n,data=data,reindex=FALSE,...)
  if (!is.null(nrow(val))) {
    val <- colSums(val)
  }
  val
}

###}}} gaussian

###{{{ gaussian variants

## Maximum Likelihood with numerical gradient + hessian

gaussian0_objective.lvm <- gaussian_objective.lvm

gaussian2_method.lvm <- "NR"
gaussian2_objective.lvm <- gaussian_objective.lvm
gaussian2_gradient.lvm <- gaussian_gradient.lvm
gaussian2_hessian.lvm <- function(p,model,n,data,...) {
##   if (FALSE) {
##   mp <- modelVar(x,p)       
##   C <- mp$C
##   iC <- Inverse(C)
##   xi <- mp$xi
##   D <- deriv(x, p=p, mom=mp,meanpar=TRUE,mu=NULL)
##   score0 <- -1/2*as.vector(iC)%*%D$dS
##   for (i in 1:nrow(data)) {
##     z <- as.numeric(mydata[i,])
##     u <- z-xi
##     score <- rbind(score,
##                    as.numeric(score0 + crossprod(u,iC)%*%D$dxi +
##                               1/2*as.vector(iC%*%tcrossprod(u)%*%iC)%*%D$dS))
##   } 
##   S <- colSums(score)
## }
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
#  val <- score(model,p=p,n=n,data=data,...)
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

