###{{{ nlminb

nlminb2 <- function(start,objective,gradient,hessian,...) {
  nlminbcontrols <- c("eval.max","iter.max","trace","abs.tol","rel.tol","x.tol","step.min")
  dots <- list(...)
  control <- list(...)$control
  control <- control[names(control)%in%nlminbcontrols]
  dots$control <- control
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

###{{{ Newton-Raphson/Scoring

NR <- function(start,objective,gradient,hessian,debug=FALSE,...) {
  dots <- list(...)
  control <- dots$control
  trace <- control$trace

##  print(control)
  if (trace>0)
  cat("Iter=0;\t\n",
      "\tp=", paste(formatC(start), collapse=" "),"\n",sep="")
  
  
  oneiter <- function(p.orig) {
    I <- hessian(p.orig)
    D <- attributes(I)$grad
    if (is.null(D)) {      
      D <- gradient(p.orig)
    } else {
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
    return(list(p=p.orig - control$gamma*iI%*%D,D=D))
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
    thetacur <- newpar$p
    if (!is.null(control$ngamma)) {
      if (control$ngamma<=gammacount) {
        control$gamma <- sqrt(control$gamma)
        gammacount <- 0
      }
    }
    if (count2==trace) {
##      cat("control$gamma=", control$gamma, "\n")
      cat("Iter=",count, ";\tD=", paste(formatC(newpar$D), collapse=" "),"\n",sep="")
      cat("\tp=", paste(formatC(thetacur), collapse=" "),"\n",sep="")      
      count2 <- 0
    }
    if (mean(newpar$D^2)<control$S.tol) break;
##    if (frobnorm(oldpar-thetacur)<control$abs.tol) break;
  }
  res <- list(par=thetacur, iterations=count, method="NR")
  return(res)
}

###}}} Newton Raphson/Scoring
