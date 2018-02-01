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

estfun <- function(start,objective,gradient,hessian,NR=FALSE,...) {
  myobj <- function(x,...) {
    S <- gradient(x,...)
    crossprod(S)[1]
  }
  if (!missing(hessian) && !is.null(hessian)) {
    mygrad <- function(x) {
      H <- hessian(x)
      S <- gradient(x)
      2*S%*%H
    }
  } else {
    hessian <- function(x) numDeriv::jacobian(gradient,x,method=lava.options()$Dmethod)
    mygrad <- function(x) {
      H <- hessian(x)
      S <- gradient(x)
      2*S%*%H
    }
  }
  if (NR) {
    op <- lava::NR(start,gradient=gradient,hessian=hessian,...)
  } else {
    op <- nlminb2(start,myobj,mygrad,hessian=NULL,...)
  }
  return(op)
}

estfun0 <- function(...,hessian=NULL) estfun(...,hessian=hessian)

###}}}

###{{{ Newton-Raphson/Scoring

##' @export
NR <- function(start,objective=NULL,gradient,hessian=NULL,debug=FALSE,control,...) {
  control0 <- list(trace=0,
                   gamma=1,
                   lambda=0,
                   ngamma=0,
                   gamma2=0,
                   backtrack=TRUE,
                   iter.max=200,
                   tol=1e-9,
                   stabil=FALSE,
                   epsilon=1e-9)
  if (!missing(control)) {
    control0[names(control)] <- control
  }
  

  ## conditions to select the step length
  if(control0$backtrack[1] == "armijo"){
    control0$backtrack <- c(1e-4,0) # page 33
  }
  if(control0$backtrack[1] == "curvature"){
    control0$backtrack <- c(0,0.9) # page 34
  }
  if(control0$backtrack[1] == "wolfe"){
      control0$backtrack <- c(1e-4,0.9)
  }
  if(!is.logical(control0$backtrack) || length(control0$backtrack)!=1){
    if(length(control0$backtrack) != 2){
      stop("control$backtrack must have length two if not TRUE or FALSE \n")
    }
    if(any(!is.numeric(control0$backtrack)) || any(abs(control0$backtrack)>1)){
      stop("elements in control$backtrack must be in [0,1] \n")
    }
    if(control0$backtrack[2]==0){
      control0$backtrack[2] <- +Inf # no Wolfe condition
    }
  }
  
  if (control0$trace>0) {
      cat("\nIter=0")
      if (!is.null(objective))
          cat("Objective=",objective(as.double(start)))
      cat(";\t\n \tp=", paste0(formatC(start), collapse=" "),"\n")
  }
  
  gradFun = !is.null(gradient)
  if (!gradFun & is.null(hessian)) {
    hessian <- function(p) {
      ff <- objective(p)
      res <- attributes(ff)$hessian
      attributes(res)$grad <- as.vector(attributes(ff)$grad)
      return(res)
    }
  }
  oneiter <- function(p.orig,Dprev,return.mat=FALSE,iter=1) {
    if (is.null(hessian)) {
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
        I <- I+control0$gamma2*sigma*diag(nrow=nrow(I))
      } else {
        sigma <- ((D)%*%t(D))
        I <- I+control0$gamma2*(sigma)
      }
    }

    Delta = control0$gamma*solve(I, D)
    iI <- Inverse(I, symmetric=TRUE, tol=control0$epsilon)
    ## Delta = control0$gamma*iI%*%as.vector(D)    
    Lambda <- 1
    if (identical(control0$backtrack, TRUE)) {
      mD0 <- mean(Dprev^2)
      mD <- mean(D^2)
      p <- p.orig + as.vector(Lambda*Delta)
      while (mD>=mD0) {
        if (gradFun) {
          D = gradient(p)
        } else {
          DI <- oneiter(p,return.mat=TRUE)
          D = DI$D
        }
        mD = mean(D^2)
        if (is.nan(mD)) mD=mD0
        Lambda <- Lambda/2
        if (Lambda<1e-4) break;
        p <- p.orig + as.vector(Lambda*Delta)
      }
      
    } else if(identical(control0$backtrack, FALSE)) {
      p <- p.orig + Lambda*Delta
    } else {  # objective(p.orig) - objective(p) <= mu*Lambda*gradient(p.orig)*Delta
      
        ## curvature
        c_D.origin_Delta <- control0$backtrack * c(rbind(D) %*% Delta)
        objective.origin <- objective(p.orig)
        p <- p.orig + as.vector(Lambda*Delta)
           
        mD0 <- c(objective.origin + Lambda * c_D.origin_Delta[1], abs(c_D.origin_Delta[2]))#    
        mD <- c(objective(p), abs(gradient(p) %*% Delta))
        count <- 0 
        while (any(mD>mD0) || any(is.nan(mD))) {
            count <- count+1
            Lambda <- Lambda/2
            if (Lambda<1e-4) break;
            p <- p.orig + Lambda*Delta
            if(!is.infinite(mD0[1])){
                mD0[1] <- objective.origin + Lambda * c_D.origin_Delta[1]#  
                mD[1] <- objective(p)
            }
            if(!is.infinite(mD0[2])){
                mD[2] <- abs(gradient(p) %*% Delta)
            }
        }
    } 
    
    return(list(p=p,D=D,iI=iI))
  }
  
  count <- count2 <- 0
  thetacur <- start
  gammacount <- 0
  Dprev <- rep(Inf,length(start))
  for (jj in seq_len(control0$iter.max)) {
    gammacount <- gammacount+1
    count <-  count+1
    count2 <- count2+1
    oldpar <- thetacur
    newpar <- oneiter(thetacur,Dprev,iter=jj)
    Dprev <- newpar$D
    thetacur <- newpar$p
    if (!is.null(control0$ngamma) && control0$ngamma>0) {
      if (control0$ngamma<=gammacount) {
        control0$gamma <- sqrt(control0$gamma)
        gammacount <- 0
      }
    }
    if (count2==control0$trace) {
        cat("Iter=", count)
        if (!is.null(objective))
            cat("Objective=",objective(as.double(newpar$p)))
        cat(";\n\tD=", paste0(formatC(newpar$D), 
                               collapse = " "), "\n")
        cat("\tp=", paste0(formatC(thetacur), collapse = " "), 
            "\n")
        count2 <- 0
    }
    if (mean(newpar$D^2)<control0$tol) break;
  }
  res <- list(par=as.vector(thetacur), iterations=count, method="NR",
              gradient=newpar$D, iH=newpar$iI)
  return(res)
}

###}}} Newton Raphson/Scoring
