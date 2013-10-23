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

  myobj <- function(...) {
    S <- gradient(...)
    return(crossprod(S)[1])
  }
  myobj <- objective
  g <- function(...) gradient(...)
  h <- function(...) numDeriv::jacobian(gradient,...)
  
  theta <- start
  print(theta)
  count <- 0
  for (i in 1:50) {
    count <- count+1
    I <- h(theta);
    I1 <- solve(I)    
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
      H <- numDeriv::jacobian(myfun,x,method=lava.options()$Dmethod)
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
    H <- numDeriv::jacobian(myfun,x,method=lava.options()$Dmethod)
    S <- gradient(x)
    2*S%*%H    
  }
  nlminb2(start,myobj,mygrad,hessian=NULL,...)
}

###}}} 

###{{{ Newton-Raphson/Scoring

NR <- function(start,objective,gradient,hessian,debug=FALSE,control,...) {
    control0 <- list(trace=0,gamma=1,lambda=0,ngamma=0,gamma2=0,
                     iter.max=200,tol=1e-15,stabil=FALSE,epsilon=1e-15)
    if (!missing(control)) {
        control0[names(control)] <- control
    }
    control <- control0
    trace <- control$trace

    if (trace>0)
        cat("\nIter=0;\t\n",
            "\tp=", paste(formatC(start), collapse=" "),"\n",sep="")

    if (missing(gradient) & missing(hessian)) {
        hessian <- function(p) {
            ff <- objective(p)
            res <- attributes(ff)$hessian
            attributes(res)$grad <- as.vector(attributes(ff)$grad)
            return(res)
        }
    }
    oneiter <- function(p.orig) {
        if (is.null(hessian)) {
            cat(".")
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
                I <- I+control$gamma2*sigma*diag(nrow(I))
            } else {
                sigma <- ((D)%*%t(D))
                I <- I+control$gamma2*(sigma)
            }
        }
        svdI <- svd(I); svdI$d0 <- numeric(length(svdI$d));
        svdI$d0[abs(svdI$d)>control$epsilon] <-
            1/svdI$d[abs(svdI$d)>control$epsilon]
        iI <- with(svdI,  (v)%*%diag(d0,nrow=length(d0))%*%t(u))   
        return(list(p=p.orig - control$gamma*iI%*%D,D=D,iI=iI))
    } 
    
    count <- count2 <- 0  
    thetacur <- start
    gammacount <- 0
    for (jj in seq_len(control$iter.max)) {
        gammacount <- gammacount+1
        count <-  count+1
        count2 <- count2+1
        oldpar <- thetacur
        newpar <- oneiter(thetacur)
        thetacur <- newpar$p
        if (!is.null(control$ngamma) && control$ngamma>0) {
            if (control$ngamma<=gammacount) {
                control$gamma <- sqrt(control$gamma)
                gammacount <- 0
            }
        }
        if (count2==trace) {
            cat("Iter=",count, ";\n\tD=", paste(formatC(newpar$D), collapse=" "),"\n",sep="")
            cat("\tp=", paste(formatC(thetacur), collapse=" "),"\n",sep="")      
            count2 <- 0
        }
        if (mean(newpar$D^2)<control$tol) break;
    }
    res <- list(par=as.vector(thetacur), iterations=count, method="NR", gradient=newpar$D, iH=newpar$iI)
    return(res)
}

###}}} NR

###{{{ NR 2

NR <- function(start,objective,gradient,hessian,debug=FALSE,control,...) {
    control0 <- list(trace=0,gamma=1,lambda=0,ngamma=0,gamma2=0,
                     iter.max=200,tol=1e-15,stabil=FALSE,epsilon=1e-15)
    if (!missing(control)) {
        control0[names(control)] <- control
    }
    trace <- control0$trace

    if (trace>0)
        cat("\nIter=0;\t\n",
            "\tp=", paste(formatC(start), collapse=" "),"\n",sep="")

    gradFun = !missing(gradient)    
    if (!gradFun & missing(hessian)) {
        hessian <- function(p) {
            ff <- objective(p)
            res <- attributes(ff)$hessian
            attributes(res)$grad <- as.vector(attributes(ff)$grad)
            return(res)
        }
    }
    oneiter <- function(p.orig,return.mat=FALSE) {
        if (is.null(hessian)) {
            cat(".")
            I <- -numDeriv::jacobian(gradient,p.orig,method=lava.options()$Dmethod)
        } else {
            I <- -hessian(p.orig)
        }
        D <- attributes(I)$grad
        if (is.null(D)) {      
            D <- gradient(p.orig)
        }
        if (return.mat) return(list(D=D,I=I))
        if (control0$stabil) {
            if (control0$lambda!=0) {
                if (control0$lambda<0) {
                    sigma <- (t(D)%*%(D))[1]          
                } else {
                    sigma <- control0$lambda
                }
                sigma <- min(sigma,10)
                I <- I+control0$gamma2*sigma*diag(nrow(I))
            } else {
                sigma <- ((D)%*%t(D))
                I <- I+control0$gamma2*(sigma)
            }
        }
        svdI <- svd(I); svdI$d0 <- numeric(length(svdI$d));
        svdI$d0[abs(svdI$d)>control0$epsilon] <-
            1/svdI$d[abs(svdI$d)>control0$epsilon]
        iI <- with(svdI,  (v)%*%diag(d0,nrow=length(d0))%*%t(u))
        Delta = control0$gamma*iI%*%D
        mD = mD0 = mean(D^2)
        Lambda <- 1
        while (mD>=mD0) {
            p <- p.orig + Lambda*Delta
            if (gradFun) {
                D = gradient(p)                
            } else {
                DI <- oneiter(p,return.mat=TRUE)
                D = DI$D
            }
            mD = mean(D^2)
            if (is.nan(mD)) mD=mD0
            Lambda <- Lambda/2
            if (Lambda<1e-12) break;
        }        
        return(list(p=p,D=D,iI=iI))
    } 
    
    count <- count2 <- 0  
    thetacur <- start
    gammacount <- 0
    for (jj in seq_len(control0$iter.max)) {
        gammacount <- gammacount+1
        count <-  count+1
        count2 <- count2+1
        oldpar <- thetacur
        newpar <- oneiter(thetacur)
        thetacur <- newpar$p
        if (!is.null(control0$ngamma) && control0$ngamma>0) {
            if (control0$ngamma<=gammacount) {
                control0$gamma <- sqrt(control0$gamma)
                gammacount <- 0
            }
        }
        if (count2==trace) {
            cat("Iter=",count, ";\n\tD=", paste(formatC(newpar$D), collapse=" "),"\n",sep="")
            cat("\tp=", paste(formatC(thetacur), collapse=" "),"\n",sep="")      
            count2 <- 0
        }
        if (mean(newpar$D^2)<control0$tol) break;
    }
    res <- list(par=as.vector(thetacur), iterations=count, method="NR",
                gradient=newpar$D, iH=newpar$iI)
    return(res)
}


###}}} Newton Raphson/Scoring
