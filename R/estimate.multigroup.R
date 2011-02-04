###{{{ estimate.multigroup

`estimate.multigroup` <- function(x, control=list(),
                                  estimator="gaussian",
                                  weight, weightname,
                                  weight2, 
                                  cluster=NULL,
                                  silent=lava.options()$silent,
                                  quick=FALSE,
                                  ...) {
  cl <- match.call()
  optim <- list(
                iter.max=lava.options()$iter.max,
                trace=ifelse(lava.options()$debug,3,0),
                gamma=1,
                gamma2=1,
                lambda=0.05,
                abs.tol=1e-9,
                epsilon=1e-10,
                delta=1e-10,
                S.tol=1e-6,
                stabil=FALSE,
                start=NULL,
                constrain=lava.options()$constrain,
                method=NULL,
                starterfun=startvalues,
                information="E",
                meanstructure=TRUE,
                sparse=FALSE,
                lbound=1e-9,
                reindex=FALSE,
                tol=1e-9)

  defopt <- lava.options()[]
  defopt <- defopt[intersect(names(defopt),names(optim))]
  optim[names(defopt)] <- defopt
  
  if (length(control)>0) {
    optim[names(control)] <- control
  }
   
  Debug("Start values...")
  if (!is.null(optim$start) & length(optim$start)==(x$npar+x$npar.mean)) {
    mystart <- optim$start
  } else {
    if (!silent) cat("Obtaining starting value...")
    mystart <- with(optim, starter.multigroup(x,meanstructure=meanstructure,starterfun=starterfun,silent=FALSE,fix=FALSE))
    if (!is.null(optim$start)) {
      pname <- names(optim$start)
      ppos <- parpos.multigroup(x,p=pname,mean=TRUE)
      if (any(!is.na(ppos)))
        mystart[match(pname,ppos)] <- optim$start[na.omit(match(ppos,pname))]
    }
    if (!silent) cat("\n")
  }
  Debug(mystart)
  Debug("Constraints...")
  ## Setup optimization constraints
  lower <- rep(-Inf, x$npar);
  for (i in 1:x$ngroup) {
    vpos <- sapply(x$parlist[[i]][variances(x$lvm[[i]],mean=FALSE)], function(y) as.numeric(substr(y,2,nchar(y))))
    if (length(vpos)>0)
    lower[vpos] <- optim$lbound  
  }
  if (optim$meanstructure)
    lower <- c(rep(-Inf,x$npar.mean), lower)
  if (any(optim$constrain)) {
    if (length(optim$constrain)!=length(lower))
      constrained <- is.finite(lower)
    else
      constrained <- optim$constrain
    lower[] <- -Inf
    optim$constrain <- TRUE    
    mystart[constrained] <- log(mystart[constrained])
  }

  if (!missing(weight)) {
    if (is.character(weight)) {
      stweight <- weight
      weight <- list()
      for (i in 1:length(x$data)) {
        newweight <- as.matrix(x$data[[i]][,stweight])
        colnames(newweight) <- index(x$lvm[[i]])$endogenous[seq_len(ncol(newweight))]
        weight <- c(weight, list(newweight))
      }
    }
  } else {
    weight <- NULL
  }
  if (!missing(weight2)) {
    if (is.character(weight2)) {      
      stweight2 <- weight2
      weight2 <- list()
      for (i in 1:length(x$data)) {
        newweight <- as.matrix(x$data[[i]][,stweight2,drop=FALSE])
        dropcol <- apply(newweight,2,function(x) any(is.na(x)))
        newweight <- newweight[,!dropcol,drop=FALSE]
        colnames(newweight) <- index(x$lvm[[i]])$endogenous[seq_len(ncol(newweight))]
        weight2 <- c(weight2, list(newweight))
      }
    }
  } else {
    weight2 <- NULL
  }

### Run hooks (additional lava plugins)
  myhooks <- gethook()
  newweight <- list()
  newweight2 <- list()
  newoptim <- newestimator <- NULL
  for (f in myhooks) {
    for ( i in 1:x$ngroup) {
      res <- do.call(f, list(x=x$lvm[[i]],data=x$data[[i]],weight=weight[[i]],weight2=weight2[[i]],estimator=estimator,optim=optim))
      if (!is.null(res$x)) x$lvm[[i]] <- res$x
      if (!is.null(res$data)) x$data[[i]] <- res$data
      newweight <- c(newweight,list(res$weight))
      newweight2 <- c(newweight2,list(res$weight2))
      if (!is.null(res$optim)) newoptim <- res$optim
      if (!is.null(res$estimator)) newestimator <- res$estimator
    }
    if (!is.null(newestimator)) estimator <- newestimator
    if (!is.null(newoptim)) optim <- newoptim
    if (!any(unlist(lapply(newweight,is.null)))) {
      weight <- newweight
    } 
    if (!any(unlist(lapply(newweight2,is.null)))) {
      weight2 <- newweight2
    }    
  }
  Method <-  paste(estimator, "_method", ".lvm", sep="")
  if (!exists(Method))
    Method <- "nlminb1"
  else
    Method <- get(Method)
  if (is.null(optim$method)) optim$method <- Method

  ## Check for random slopes
  xXx <- exogenous(x)
  Xfix <- FALSE
  Xconstrain <- FALSE
  xfix <- list()
  xconstrain <- list
  for (i in 1:x$ngroup) {
    x0 <- x$lvm[[i]]
    data0 <- x$data[[i]]
    xfix0 <- colnames(data0)[(colnames(data0)%in%parlabels(x0,exo=TRUE))]
    xconstrain0 <- intersect(unlist(lapply(constrain(x0),function(z) attributes(z)$args)),manifest(x0))
    xconstrain <- c(xconstrain,list(xconstrain0))
    xfix <- c(xfix, list(xfix0))
    if (length(xfix0)>0) { ## Yes, random slopes
      Xfix<-TRUE
    }
    if (length(xconstrain0)>0) {
      Xconstrain <- TRUE
    }
  }
  
  ## Non-linear parameter constraints involving observed variables? (e.g. nonlinear regression)

  ## Define objective function and first and second derivatives
  ObjectiveFun  <- paste(estimator, "_objective", ".lvm", sep="")
  GradFun  <- paste(estimator, "_gradient", ".lvm", sep="")
  if (!exists(ObjectiveFun) & !exists(GradFun)) stop("Unknown estimator.")

  InformationFun <- paste(estimator, "_hessian", ".lvm", sep="")
  
  
  ##parord <- modelPar(x,1:length(mystart),debug=debug)$p
  parord <- modelPar(x,1:with(x,npar+npar.mean))$p
  mymodel <- x

  parkeep <- c()
  myclass <- c("multigroupfit","lvmfit")  
  myfix <- list()

  if (Xfix | Xconstrain) { ## Model with random slopes:
################################################################################
################################################################################
################################################################################

    if (Xfix) {
      myclass <- c(myclass,"lvmfit.randomslope")
      for (k in 1:x$ngroup) {
        x1 <- x0 <- x$lvm[[k]]
        data0 <- x$data[[k]]
        
        nrow <- length(vars(x0))
        xpos <- lapply(xfix[[k]],function(y) which(regfix(x0)$labels==y))
        colpos <- lapply(xpos, function(y) ceiling(y/nrow))
        rowpos <- lapply(xpos, function(y) (y-1)%%nrow+1)
        myfix0 <- list(var=xfix[[k]], col=colpos, row=rowpos)
        myfix <- c(myfix, list(myfix0))
      
        for (i in 1:length(myfix0$var))
          for (j in 1:length(myfix0$col[[i]])) 
            regfix(x0,
                   from=vars(x0)[myfix0$row[[i]][j]],to=vars(x0)[myfix0$col[[i]][j]]) <-
                     colMeans(data0[,myfix0$var[[i]],drop=FALSE],na.rm=TRUE)
        index(x0) <- reindex(x0,zeroones=TRUE,deriv=TRUE)
        x$lvm[[k]] <- x0
        yvars <- endogenous(x0)
        parkeep <- c(parkeep, parord[[k]][coef(x1,mean=TRUE)%in%coef(x0,mean=TRUE)])
      }
      parkeep <- sort(unique(parkeep))
      ## Alter start-values:

      if (length(mystart)!=length(parkeep))
        mystart <- mystart[parkeep]
      lower <- lower[parkeep]
      x <- multigroup(x$lvm,x$data,fix=FALSE,exo.fix=FALSE) 
    }  
    
    parord <- modelPar(x,1:length(mystart))$p    
    mydata <- list()
    for (i in 1:x$ngroup) {      
      mydata <- c(mydata, list(as.matrix(x$data[[i]][,manifest(x$lvm[[i]])])))
    }
    
    myObj <- function(theta) {
      if (optim$constrain)
        theta[constrained] <- exp(theta[constrained])
      pp <- modelPar(x,theta)$p
      res <- 0
    for (k in 1:x$ngroup) {
##      for (k in 5:5) {
        x0 <- x$lvm[[k]]
        data0 <- x$data[[k]]
        if (Xfix) {
          xfix0 <- xfix[[k]]
          myfix0 <- myfix[[k]]
        }
        p0 <- pp[[k]]
        myfun <- function(ii) {
          if (Xfix) 
          for (i in 1:length(myfix0$var)) {
            x0$fix[cbind(myfix0$row[[i]],myfix0$col[[i]])] <-
              index(x0)$A[cbind(myfix0$row[[i]],myfix0$col[[i]])] <-
                data0[ii,xfix0[i]]
          }
          if (is.list(weight2[[k]][ii,])) {
            res <- do.call(ObjectiveFun, list(x=x0, p=p0,
                                              data=data0[ii,manifest(x0),drop=FALSE],
                                              n=1, S=NULL, weight=weight[[k]][ii,],
                                              weight2=weight2[[k]]))
            
          } else {
            res <- do.call(ObjectiveFun, list(x=x0, p=p0,
                                              data=data0[ii,manifest(x0),drop=FALSE],
                                              n=1, S=NULL, weight=weight[[k]][ii,],
                                              weight2=weight2[[k]][ii,]))
          }
          return(res)
        }
        res <- res + sum(sapply(1:nrow(mydata[[k]]),myfun))
      }
      res
    }

    myGrad <- function(theta) {
      if (optim$constrain) {
        theta[constrained] <- exp(theta[constrained])
      }
      pp <- modelPar(x,theta)$p
      D0 <- res <- rbind(numeric(length(mystart)))
      for (k in 1:x$ngroup) {
        if (Xfix) {
          myfix0 <- myfix[[k]]
        }
        x0 <- x$lvm[[k]]
        myfun <- function(ii) {
          if (Xfix)
          for (i in 1:length(myfix0$var)) {
            x0$fix[cbind(myfix0$row[[i]],myfix0$col[[i]])] <-
              index(x0)$A[cbind(myfix0$row[[i]],myfix0$col[[i]])] <-
                x$data[[k]][ii,xfix[[k]][i]]
          }
          if (is.list(weight2[[k]][ii,])) {
            
          } else {
            val <- do.call(GradFun, list(x=x0, p=pp[[k]],
                                         data=mydata[[k]][ii,,drop=FALSE], n=1,
                                         S=NULL,
                                         weight=weight[[k]][ii,],
                                         weight2=weight2[[k]][ii,]))
          }
          return(val)
        }
        D <- D0; D[parord[[k]]] <- rowSums(sapply(1:nrow(mydata[[k]]),myfun))
        res <- res+D
      }
      if (optim$constrain) {
        res[constrained] <- res[constrained]*theta[constrained]
      }
      return(as.vector(res))
    }

    myInformation <- function(theta) {
      theta0 <- theta
      if (optim$constrain) {
        theta[constrained] <- exp(theta[constrained])
      }
      pp <- modelPar(x,theta)$p
      I0 <- res <- matrix(0,length(theta),length(theta))
      grad <- grad0 <- numeric(length(theta))
      for (k in 1:x$ngroup) {
        x0 <- x$lvm[[k]]
        if (Xfix) {
          myfix0 <- myfix[[k]]
        }
        myfun <- function(ii) {
          if (Xfix)
          for (i in 1:length(myfix0$var)) {
            x0$fix[cbind(myfix0$row[[i]],myfix0$col[[i]])] <- index(x0)$A[cbind(myfix0$row[[i]],myfix0$col[[i]])] <-
              x$data[[k]][ii,xfix[[k]][i]]
          }
          I <- I0
          J <- do.call(InformationFun,
                       list(x=x0, p=pp[[k]],
                            data=mydata[[k]][ii,], n=1,
                            S=NULL,
                            weight=weight[[k]][ii,],
                            weight2=weight2[[k]][ii,],
                            type=optim$information
                            )
                       )
          D <- grad0
          if (!is.null(attributes(J)$grad)) {
            D[ parord[[k]] ] <- attributes(J)$grad
            attributes(I)$grad <- D
          }
          I[ parord[[k]], parord[[k]] ] <- J
##                                                   list(p=pp[[k]], model=x0, data=x$data[[k]], n=1, weight=weight[[k]][ii,]))
          return(I)
        }      
        L <- lapply(1:nrow(x$data[[k]]),function(x) myfun(x))
        if (!is.null(attributes(L[[1]])$grad))
          grad <- grad + rowSums(matrix((unlist(lapply(L,function(x) attributes(x)$grad))),ncol=length(L)))        
        res <- res + apply(array(unlist(L),dim=c(length(theta),length(theta),nrow(x$data[[k]]))),c(1,2),sum)
      }
      if (!is.null(attributes(L[[1]])$grad))
        attributes(res)$grad <- grad
      return(res)   
    }
  } else { ## Model without random slopes:    
################################################################################
#########################################2#######################################
################################################################################
    
    myObj <- function(theta) {     
      theta0 <- theta
      if (optim$constrain) {
        theta[constrained] <- exp(theta[constrained])
      }
      pp <- modelPar(x,theta)$p
      res <- c()
      for (i in 1:x$ngroup) { 
        res <- c(res,
                 with(x$samplestat[[i]],do.call(ObjectiveFun, list(x=x$lvm[[i]], p=pp[[i]], data=x$data[[i]][,index(x$lvm[[i]])$manifest,drop=FALSE], S=S, mu=mu, n=n, weight=weight[[i]], weight2=weight2[[i]])))
                 )
      }
        sum(res)
    }
    
    if (!exists(GradFun)) {
      myGrad <- NULL
    } else  {
      myGrad <- function(theta) {
        theta0 <- theta
        if (optim$constrain) {
          theta[constrained] <- exp(theta[constrained])
        }
        pp <- modelPar(x,theta)$p
        D0 <- res <- rbind(numeric(length(theta)))
        for (i in 1:x$ngroup) {
          repval <- with(x$samplestat[[i]],
                         do.call(GradFun, list(x=x$lvm[[i]],p=pp[[i]],
                                               data=x$data[[i]][,index(x$lvm[[i]])$manifest,drop=FALSE],
                                               S=S,mu=mu,n=n,
                                               weight=weight[[i]], weight2=weight2[[i]])))
          D <- D0; D[ parord[[i]] ] <- repval
        res <- res + D
        }
        if (optim$constrain) {
          res[constrained] <- res[constrained]*theta[constrained]
        }
        return(as.vector(res))
      }
    }
    
    myInformation <- function(theta) {
      theta0 <- theta
      if (optim$constrain) {
        theta[constrained] <- exp(theta[constrained])
      }
      ##      parord <- modelPar(x,1:length(theta),debug=debug)$p
      pp <- modelPar(x,theta)$p
      I0 <- res <- matrix(0,length(theta),length(theta))
      for (i in 1:x$ngroup) {
        I <- I0;
        I[ parord[[i]], parord[[i]] ] <- with(x$samplestat[[i]], do.call(InformationFun, list(p=pp[[i]], x=x$lvm[[i]], data=x$data[[i]],
                                                                                              S=S, mu=mu, n=n, weight=weight[[i]],
                                                                                              weight2=weight2[[i]],
                                                                                              type=optim$information)))
        ##with(x$samplestat[[i]], information(x$lvm[[i]],p=pp[[i]],n=n))
        res <- res + I
      }
      D <- myGrad(theta0)
      if (optim$constrain) {
        res[constrained,-constrained] <- apply(res[constrained,-constrained,drop=FALSE],2,function(x) x*theta[constrained]);
        res[-constrained,constrained] <- t(res[constrained,-constrained])
        if (sum(constrained)==1) {
          res[constrained,constrained] <- res[constrained,constrained]*outer(theta[constrained],theta[constrained]) - (D[constrained])
        } else {
          res[constrained,constrained] <- res[constrained,constrained]*outer(theta[constrained],theta[constrained]) - diag(D[constrained])
        }
      }
      attributes(res)$grad <- D
      return(res)
    }    
  }
  
################################################################################
################################################################################
################################################################################

  if (!exists(InformationFun)) myInformation <- NULL
  else if (is.null(get(InformationFun))) myInformation <- NULL
  if (is.null(get(GradFun))) myGrad <- NULL
    
  if (!silent) cat("Optimizing objective function...")
  if (lava.options()$debug) {
    print(lower)
    print(optim$constrain)
    print(optim$method)
  }

  opt <- do.call(optim$method,
                 list(start=mystart, objective=myObj, gradient=myGrad, hessian=myInformation, lower=lower, control=optim))
  if (!silent) cat("\n")

  opt$estimate <- opt$par
  if (optim$constrain) {
    opt$estimate[constrained] <- exp(opt$estimate[constrained])
  }  
  if (quick) return(list(opt=opt,vcov=NA))
  
  if (is.null(myGrad)) {
    if (!require("numDeriv")) {
      opt$gradient <- naiveGrad(myObj, opt$estimate)
    } else {
      opt$gradient <- grad(myObj, opt$estimate, method="Richardson")
    }
  } else {
    opt$gradient <- myGrad(opt$estimate)
  }

  if (is.null(myInformation)) {
    if (!require("numDeriv")) stop("I do not know how to calculate the asymptotic variance of this estimator.
For numerical approximation please install the library 'numDeriv'.")
    ##    cat("Using a numerical approximation of hessian...\n");
    if (!is.null(myGrad))
      myInformation <- function(theta) jacobian(myGrad, theta, method=lava.options()$Dmethod)
    else
      myInformation <- function(theta) -hessian(myObj, theta, method=lava.options()$Dmethod)
  }
  I <- myInformation(opt$estimate)
  asVar <- tryCatch(solve(I),
                    error=function(e) matrix(NA, length(mystart), length(mystart)))

##  if (!silent) cat("\n")

    
  res <- list(model=x, model0=mymodel, call=cl, opt=opt, meanstructure=optim$meanstructure, vcov=asVar, estimator=estimator, weight=weight, weight2=weight2, cluster=cluster)
  class(res) <- myclass

  myhooks <- gethook("post.hooks")
  for (f in myhooks) {
    res0 <- do.call(f,list(x=res))
    if (!is.null(res0))
      res <- res0
  }

  return(res)
}

###}}}

###{{{ estimate.list

`estimate.list` <-
function(x, data, silent=lava.options()$silent, fix, missing=FALSE,  ...) {
  if (missing(data)) {
    return(estimate(x[[1]],x[[2]],missing=missing,...))
  }  
  nm <- length(x)
  if (nm==1) {
    return(estimate(x[[1]],data,missing=missing,...))
  }
  if (!all(unlist(lapply(x, function(y) class(y)[1]=="lvm")))) stop ("Expected a list of 'lvm' objects.")
  if (is.data.frame(data)) {
    warning("Only one dataset - going for standard analysis on each submodel.")
    res <- c()
    for (i in 1:nm) {
      res <- c(res, list(estimate(x[[i]],data=data,silent=TRUE,missing=missing, ...)))
    }
    return(res)
  }

  if (nm!=length(data)) stop("Supply dataset for each model")

  Xfix <- FALSE
  xfix <- list()
  for (i in 1:length(x)) {
    data0 <- data[[i]]
    xfix0 <- colnames(data0)[(colnames(data0)%in%parlabels(x[[i]],exo=TRUE))]
    xfix <- c(xfix, list(xfix0))
    if (length(xfix0)>0) { ## Yes, random slopes
      Xfix<-TRUE
    }
  }
  if (missing(fix)) {
    fix <- ifelse(Xfix,FALSE,TRUE)
  }  
  
  mg <- multigroup(x,data,fix=fix,missing=missing,...)
  res <- estimate(mg,...)
    
  return(res) 
}

###}}}z
