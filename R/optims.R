###{{{ nlminb

nlminb2 <- function(start,objective,gradient,hessian,...) {
  nlminbcontrols <- c("eval.max","iter.max","trace","abs.tol","rel.tol","x.tol","step.min")
  dots <- list(...)
  control <- list(...)$control
  control <- control[names(control)%in%nlminbcontrols]
  dots$control <- control
  if (length(dots$trace)>0 && dots$trace>0) cat("\n")
  mypar <- c(list(start=start,objective=objective,gradient=gradient,hessian=hessian),dots)
  mypar["debug"] <- NULL
  do.call("nlminb", mypar)
##  nlminb(start,objective,gradient=NULL,hessian=NULL,control=control,...)
}

nlminb1 <- function(start,objective,gradient,hessian,...) {
  nlminb2(start,objective,gradient=gradient,hessian=NULL,...)
}
nlminb0 <- function(start,objective,gradient,hessian,...) {
  nlminb2(start,objective,gradient=NULL,hessian=NULL,...)
}

###}}} nlminb

###{{{ estfun


### just for testing
estfun2 <- function(start,objective,gradient,hessian,...) {

  browser()
  myobj <- function(...) {
##    return(objective(...))
##    return(objective(...))
    S <- gradient(...)
    return(crossprod(S)[1])
  }
  myobj <- objective
  ##g <- function(...) numDeriv::grad(myobj,...)
  g <- function(...) gradient(...)
##  h <- function(...) numDeriv::hessian(myobj,...)
  h <- function(...) numDeriv::jacobian(gradient,...)
  ##nlminb(start,myobj,g,h,control=list(trace=1))
  
  theta <- start
  print(theta)
  count <- 0
  for (i in 1:50) {
    count <- count+1
    I <- h(theta);
    I1 <- solve(I)    
    ##E <- eigen(I); L <- E$values;    
    ##L[L<1e-2] <- 0; L[L>0] <- 1/L[L>0]    
    ##I1 <- with(E, vectors%*%diag(L)%*%t(vectors))
    S <- g(theta)
    theta <- theta + (-0.5*I1%*%S)
    print(S)
    print(theta[,1])
  }
  
  return(theta)
    res <- list(par=theta[,1], iterations=count, method="estfun0", gradient=g(start))
 
}


estfun <- function(start,objective,gradient,hessian,...) {
  myobj <- function(...) {
    S <- gradient(...)
    crossprod(S)[1]
  }
  if (!is.null(hessian)) {
    mygrad <- function(x) {
      H <- hessian(x)
      S <- gradient(x)
      2*S%*%H    
    }    
  } else {
    mygrad <- function(x) {
      myfun <- function(z) gradient(z)
      H <- jacobian(myfun,x,method=lava.options()$Dmethod)
      S <- gradient(x)
      2*S%*%H    
    }
  }
  nlminb2(start,myobj,mygrad,hessian=NULL,...)
}

estfun0 <- function(start,objective,gradient,hessian,...) {
  myobj <- function(...) {
    S <- gradient(...)
    crossprod(S)[1]
  }
  mygrad <- function(x) {
    myfun <- function(z) gradient(z)
    H <- jacobian(myfun,x,method=lava.options()$Dmethod)
    S <- gradient(x)
    2*S%*%H    
  }
  nlminb2(start,myobj,mygrad,hessian=NULL,...)
}



###}}} 


###{{{ Newton-Raphson/Scoring

NR <- function(start,objective,gradient,hessian,debug=FALSE,...) {
  dots <- list(...)
  control <- dots$control
  trace <- control$trace

##  print(control)
  if (trace>0)
  cat("\nIter=0;\t\n",
      "\tp=", paste(formatC(start), collapse=" "),"\n",sep="")
  

  oneiter <- function(p.orig) {
    if (is.null(hessian)) {
      cat(".")
      ##      I <- numDeriv::jacobian(gradient,p.orig)      
      ##      I <- numDeriv::jacobian(gradient,p.orig,method="simple")
      I <- numDeriv::jacobian(gradient,p.orig,method=lava.options()$Dmethod)
    } else {
      I <- hessian(p.orig)
    }
    D <- attributes(I)$grad
    if (is.null(D)) {      
      D <- gradient(p.orig)
    }
    
    if (control$stabil) {
      if (control$lambda!=0) {
        if (control$lambda<0) {
          sigma <- (t(D)%*%(D))[1]          
        } else {
          sigma <- control$lambda
        }
        sigma <- min(sigma,10)
        I <- I+sigma*diag(nrow(I))
      } else {
        sigma <- ((D)%*%t(D))
        ##        K <- max(diag(sigma))
        I <- I+control$gamma2*(sigma)
      }
    }
    svdI <- svd(I); svdI$d0 <- numeric(length(svdI$d));
    ##  delta <- 0
    ##svdI$d0 <- 1/(abs(svdI$d)+delta)
    svdI$d0[abs(svdI$d)>control$epsilon] <-
      1/svdI$d[abs(svdI$d)>control$epsilon]
    ##    svdI$d0 <- 1/svdI$d0
    ##+control$delta)
    ##    Debug(list("d0",svdI$d0), debug)    
    ##    Debug(list("v",svdI$v), debug)    
    ##    save(svdI, file="I.rda")
    iI <- with(svdI,  (v)%*%diag(d0,nrow=length(d0))%*%t(u))
    ##    iI <- with(svdI,  (v)%*%diag(1/d)%*%t(u))
    ##    iI <- solve(I)
    ##    I <- I + 0.001*diag(nrow(I))
    ##browser()
    return(list(p=p.orig - control$gamma*iI%*%D,D=D,iI=iI))
  } 

  
  count <- count2 <- 0  
  thetacur <- start
  gammacount <- 0
  for (jj in 1:control$iter.max) {
    gammacount <- gammacount+1
    count <-  count+1
    count2 <- count2+1
    oldpar <- thetacur
    newpar <- oneiter(thetacur)
    ##browser()
    thetacur <- newpar$p
    if (!is.null(control$ngamma)) {
      if (control$ngamma<=gammacount) {
        control$gamma <- sqrt(control$gamma)
        gammacount <- 0
      }
    }
    if (count2==trace) {
##      cat("control$gamma=", control$gamma, "\n")
      cat("Iter=",count, ";\n\tD=", paste(formatC(newpar$D), collapse=" "),"\n",sep="")
      cat("\tp=", paste(formatC(thetacur), collapse=" "),"\n",sep="")      
      count2 <- 0
    }
    if (mean(newpar$D^2)<control$tol) break;
##    if (frobnorm(oldpar-thetacur)<control$abs.tol) break;
  }
  res <- list(par=as.vector(thetacur), iterations=count, method="NR", gradient=newpar$D, iH=newpar$iI)
  return(res)
}

###}}} Newton Raphson/Scoring
