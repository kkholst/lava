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


## Newton-Raphson/Scoring

##' @title Newton-Raphson method
##' 
##' @param start Starting value
##' @param objective Optional objective function (used for selecting step length)
##' @param gradient gradient
##' @param hessian hessian (if NULL a numerical derivative is used)
##' @param control optimization arguments (see details)
##' @param args Optional list of arguments parsed to objective, gradient and hessian
##' @param ... additional arguments parsed to lower level functions
##' @details
##' \code{control} should be a list with one or more of the following components:
##' \itemize{
##' \item{trace} integer for which output is printed each 'trace'th iteration
##' \item{iter.max} number of iterations
##' \item{stepsize}: Step size (default 1)
##' \item{nstepsize}: Increase stepsize every nstepsize iteration (from stepsize to 1)
##' \item{tol}: Convergence criterion (gradient)
##' \item{epsilon}: threshold used in pseudo-inverse
##' \item{backtrack}: In each iteration reduce stepsize unless solution is improved according to criterion (gradient, armijo, curvature, wolfe)
##' }
##' @export
##' @examples
##' # Objective function with gradient and hessian as attributes
##' f <- function(z) {
##'     x <- z[1]; y <- z[2]
##'     val <- x^2 + x*y^2 + x + y
##'     structure(val, gradient=function(x) c(2*x+y^2+1, x+1),
##'               hessian=function(x) c(2, 0))
##' }
##' NR(c(0,0),f)
##' 
##' 
##' # Parsing arguments to the function and
##' g <- function(x,y) (x*y+1)^2
##' NR(0, gradient=g, args=list(y=2), control=list(trace=1,tol=1e-20))
##' 
##' 
NR <- function(start,objective=NULL,gradient=NULL,hessian=NULL,control,args=NULL,...) {
  control0 <- list(trace=0,
                   stepsize=1,
                   lambda=0,
                   ngamma=0,
                   gamma2=0,
                   backtrack=TRUE,
                   iter.max=200,
                   tol=1e-6,
                   stabil=FALSE,
                   epsilon=1e-9)
  if (!missing(control)) {
      control0[names(control)] <- control
      # Backward compatibility:
      if (!is.null(control0$gammma)) control0$stepsize <- control0$gamma
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
  obj <- objective
  grad <- gradient
  hess <- hessian
  if (!is.null(args)) {
      if (!is.list(args)) args <- list(args)
      if (!is.null(objective))
          obj <- function(p) do.call(objective, c(list(p),args))
      if (!is.null(gradient))
          grad <- function(p) do.call(gradient, c(list(p),args))
      if (!is.null(hessian))
          hess <- function(p) do.call(hessian, c(list(p),args))
  }

  if (control0$trace>0) {
      cat("\nIter=0")
      if (!is.null(obj))
          cat("Objective=",obj(as.double(start)))
      cat(";\t\n \tp=", paste0(formatC(start), collapse=" "),"\n")
  }

  gradFun = !is.null(grad)
  if (!gradFun & is.null(hess)) {
    hess <- function(p) {
      ff <- obj(p)
      res <- attributes(ff)$hessian
      attributes(res)$grad <- as.vector(attributes(ff)$grad)
      return(res)
    }
    grad <- function(p) numDeriv::jacobian(obj,p)
    hess <- NULL
  }
  oneiter <- function(p.orig,Dprev,return.mat=FALSE,iter=1) {
    if (is.null(hess)) {
      I <- -numDeriv::jacobian(grad,p.orig,method=lava.options()$Dmethod)
    } else {
      I <- -hess(p.orig)
    }
    D <- attributes(I)$grad
    if (is.null(D)) {
      D <- grad(p.orig)
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
    iI <- Inverse(I, symmetric=TRUE, tol=control0$epsilon)
    Delta <- control0$stepsize*tryCatch(solve(I, cbind(as.vector(D))),
                            error=function(...) { ## Fall back to Pseudo-Inverse using SVD:
                                iI%*%cbind(as.vector(D))})
    Lambda <- 1 ## Initial step-size
    if (identical(control0$backtrack, TRUE)) {
      mD0 <- mean(Dprev^2)
      mD <- mean(D^2)
      p <- p.orig + as.vector(Lambda*Delta)
      while (mD>=mD0) {
        if (gradFun) {
          D = grad(p)
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
    } else {  # objective(p.orig) - obj(p) <= mu*Lambda*grad(p.orig)*Delta

        ## curvature
        c_D.origin_Delta <- control0$backtrack * c(rbind(D) %*% Delta)
        objective.origin <- obj(p.orig)
        p <- p.orig + as.vector(Lambda*Delta)

        mD0 <- c(objective.origin + Lambda * c_D.origin_Delta[1], abs(c_D.origin_Delta[2]))#
        mD <- c(obj(p), abs(grad(p) %*% Delta))
        count <- 0
        while (any(mD>mD0) || any(is.nan(mD))) {
            count <- count+1
            Lambda <- Lambda/2
            if (Lambda<1e-4) break;
            p <- p.orig + Lambda*Delta
            if(!is.infinite(mD0[1])){
                mD0[1] <- objective.origin + Lambda * c_D.origin_Delta[1]#
                mD[1] <- obj(p)
            }
            if(!is.infinite(mD0[2])){
                mD[2] <- abs(grad(p) %*% Delta)
            }
        }
    }

    return(list(p=p,D=D,iI=iI))
  }

  count <- count2 <- 0
  thetacur <- start
  stepsizecount <- 0
  Dprev <- rep(Inf,length(start))
  for (jj in seq_len(control0$iter.max)) {
    stepsizecount <- stepsizecount+1
    count <-  count+1
    count2 <- count2+1
    newpar <- oneiter(thetacur,Dprev,iter=jj)
    Dprev <- newpar$D
    thetacur <- newpar$p
    if (!is.null(control0$nstepsize) && control0$nstepsize>0) {
      if (control0$nstepsize<=stepsizecount) {
        control0$stepsize <- sqrt(control0$stepsize)
        stepsizecount <- 0
      }
    }
    if (count2==control0$trace) {
        cat("Iter=", count)
        if (!is.null(obj))
            cat("Objective=",obj(as.double(newpar$p)))
        cat(";\n\tD=", paste0(formatC(newpar$D),
                               collapse = " "), "\n")
        cat("\tp=", paste0(formatC(thetacur), collapse = " "),
            "\n")
        count2 <- 0
    }
    if (mean(newpar$D^2)^.5<control0$tol) break;
  }
  res <- list(par=as.vector(thetacur), iterations=count, method="NR",
              gradient=newpar$D, iH=newpar$iI)
  return(res)
}

