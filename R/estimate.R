##' @export
estimate <- function(x,...) UseMethod("estimate")

###{{{ estimate.lvm

##' Estimation of parameters in a Latent Variable Model (lvm)
##' 
##' Estimate parameters. MLE, IV or user-defined estimator.
##' 
##' A list of parameters controlling the estimation and optimization procedures
##' is parsed via the \code{control} argument. By default Maximum Likelihood is
##' used assuming multivariate normal distributed measurement errors. A list
##' with one or more of the following elements is expected:
##'
##' \describe{
##' \item{start:}{Starting value. The order of the parameters can be shown by
##' calling \code{coef} (with \code{mean=TRUE}) on the \code{lvm}-object or with
##' \code{plot(..., labels=TRUE)}. Note that this requires a check that it is
##' actual the model being estimated, as \code{estimate} might add additional
##' restriction to the model, e.g. through the \code{fix} and \code{exo.fix}
##' arguments. The \code{lvm}-object of a fitted model can be extracted with the
##' \code{Model}-function.}
##'
##' \item{starterfun:}{Starter-function with syntax
##' \code{function(lvm, S, mu)}.  Three builtin functions are available:
##' \code{startvalues}, \code{startvalues2}, \ codestartvalues3.}
##'
##' \item{estimator:}{ String defining which estimator to use (Defaults to
##' ``\code{gaussian}'')}
##'
##' \item{meanstructure}{Logical variable indicating
##' whether to fit model with meanstructure.}
##'
##' \item{method:}{ String pointing to
##' alternative optimizer (e.g. \code{optim} to use simulated annealing).}
##'
##' \item{control:}{ Parameters passed to the optimizer (default
##' \code{stats::nlminb}).}
##' 
##' \item{tol:}{ Tolerance of optimization constraints on lower limit of
##' variance parameters.  } }
##'
##' @aliases estimate.list estimate
##' @param x \code{lvm}-object
##' @param data \code{data.frame}
##' @param estimator String defining the estimator (see details below)
##' @param control control/optimization parameters (see details below)
##' @param weight Optional weights to used by the chosen estimator.
##' @param weightname Weight names (variable names of the model) in case
##' \code{weight} was given as a vector of column names of \code{data}
##' @param weight2 Optional second set of weights to used by the chosen
##' estimator.
##' @param cluster Vector (or name of column in \code{data}) that identifies
##' correlated groups of observations in the data leading to variance estimates
##' based on a sandwich estimator
##' @param missing Logical variable indiciating how to treat missing data.
##' Setting to FALSE leads to complete case analysis. In the other case
##' likelihood based inference is obtained by integrating out the missing data
##' under assumption the assumption that data is missing at random (MAR).
##' @param index For internal use only
##' @param graph For internal use only
##' @param fix Logical variable indicating whether parameter restriction
##' automatically should be imposed (e.g. intercepts of latent variables set to
##' 0 and at least one regression parameter of each measurement model fixed to
##' ensure identifiability.)
##' @param quick If TRUE the parameter estimates are calculated but all
##' additional information such as standard errors are skipped
##' @param silent Logical argument indicating whether information should be
##' printed during estimation
##' @param \dots Additional arguments to be passed to the low level functions
##' @return A \code{lvmfit}-object.
##' @author Klaus K. Holst
##' @seealso \code{\link{score}}, \code{\link{information}}, ...
##' @keywords models regression
##' @S3method estimate lvm
##' @method estimate lvm
##' @examples
##' 
##' m <- lvm(list(y~v1+v2+v3+v4,c(v1,v2,v3,v4)~x))
##' \donttest{
##' plot(m)
##' }
##' dd <- sim(m,10000) ## Simulate 10000 observations from model
##' e <- estimate(m, dd) ## Estimate parameters
##' e
##'
`estimate.lvm` <-
function(x, data,
         estimator="gaussian",
         control=list(),
         missing=FALSE,
         weight,
         weightname,
         weight2,
         cluster,
         fix,
         index=TRUE,
         graph=FALSE,
         silent=lava.options()$silent,
         quick=FALSE,
         ...) { 

  if (length(exogenous(x)>0)) {
    catx <- categorical2dummy(x,data)
    x <- catx$x; data <- catx$data
  }
  
  cl <- match.call()

  optim <- list(
                iter.max=lava.options()$iter.max,
                trace=ifelse(lava.options()$debug,3,0),
                gamma=1,
                gamma2=1,
                ngamma=NULL,
                lambda=0.05,
                abs.tol=1e-9,
                epsilon=1e-10,
                delta=1e-10,
                rel.tol=1e-10,
                S.tol=1e-5,
                stabil=FALSE,
                start=NULL,
                constrain=lava.options()$constrain,
                method=NULL,
                starterfun="startvalues",
                information="E",
                meanstructure=TRUE,
                sparse=FALSE,
                tol=1e-9)

  defopt <- lava.options()[]
  defopt <- defopt[intersect(names(defopt),names(optim))]
  optim[names(defopt)] <- defopt
  if (length(control)>0) {
    optim[names(control)] <- control
  }

  if (!lava.options()$exogenous) exogenous(x) <- NULL
  ## Random-slopes:
  redvar <- intersect(intersect(parlabels(x),latent(x)),colnames(data))
  if (length(redvar)>0)
    warning(paste("Remove latent variable colnames from dataset",redvar))
  xfix <- setdiff(colnames(data)[(colnames(data)%in%parlabels(x,exo=TRUE))],latent(x))
  if (missing(fix)) {
    fix <- ifelse(length(xfix)>0,FALSE,TRUE)
  }
  Debug(list("start=",optim$start))
  ## Weights...
  if (!missing(weight)) {
    if (is.character(weight)) {
      weight <- data[,weight,drop=FALSE]
      if (!missing(weightname)) {
        colnames(weight) <- weightname
      } else {
        yvar <- index(x)$endogenous
        nw <- seq_len(min(length(yvar),ncol(weight)))
        colnames(weight)[nw] <- yvar[nw]
      }      
    }
    weight <- cbind(weight)
    
  } else {
    weight <- NULL
  }
  if (!missing(weight2)) {
    if (is.character(weight2)) {
      weight2 <- data[,weight2]
    }
  } else {
    weight2 <- NULL
  }
  ## Correlated clusters...
  if (!missing(cluster)) {
    if (is.character(cluster)) {
      cluster <- data[,cluster]
    } 
  } else {
    cluster <- NULL
  }
  
  Debug("procdata")
  if (missing) { ## Remove rows with missing covariates
    xmis <- apply(data[,exogenous(x),drop=FALSE],1,function(x) any(is.na(x)))    
    data <- data[which(!xmis),]
  }
  ## e1<- estimate(m1,testdata[which(!xmis),-1],missing=TRUE)
  if (!missing & (is.matrix(data) | is.data.frame(data))) {    
    data <- na.omit(data[,intersect(colnames(data),c(manifest(x),xfix)),drop=FALSE])
  }

  dd <- procdata.lvm(x,data=data)
  S <- dd$S; mu <- dd$mu; n <- dd$n
  Debug(list("n=",n))  
  Debug(list("S=",S))
  Debug(list("mu=",mu))

  ##  if (fix)
  {
    var.missing <- setdiff(vars(x),colnames(S))
    if (length(var.missing)>0) {## Convert to latent:
      new.lat <- setdiff(var.missing,latent(x))
      if (length(new.lat)>0)
        x <- latent(x, new.lat)
    }
  }

  ## Run hooks (additional lava plugins)
  myhooks <- gethook()
  for (f in myhooks) {    
    res <- do.call(f, list(x=x,data=data,weight=weight,weight2=weight2,estimator=estimator,optim=optim))
    if (!is.null(res$x)) x <- res$x
    if (!is.null(res$data)) data <- res$data
    if (!is.null(res$weight)) weight <- res$weight
    if (!is.null(res$weight2)) weight2 <- res$weight2
    if (!is.null(res$optim)) optim <- res$optim
    if (!is.null(res$estimator)) estimator <- res$estimator
  }
 
  Method <-  paste(estimator, "_method", ".lvm", sep="")
  if (!exists(Method)) {
    Method <- "nlminb1"
  } else {
    Method <- get(Method)
  }
  if (is.null(optim$method)) optim$method <- Method

  if (!quick & index) {
    ## Proces data and setup some matrices
    x <- fixsome(x, measurement.fix=fix, S=S, mu=mu, n=n,debug=!silent)
    if (!silent)
      message("Reindexing model...\n")
    if (length(xfix)>0) {
      index(x) <- reindex(x,sparse=optim$sparse,zeroones=TRUE,deriv=TRUE)
    } else {
      x <- updatelvm(x,sparse=optim$sparse,zeroones=TRUE,deriv=TRUE,mean=TRUE)
    }
  }

  if (is.null(estimator) || estimator==FALSE) {
    return(x)
  }  
  k <- length(manifest(x))
  Debug(list("S=",S))
  if (!optim$meanstructure) {
    mu <- NULL
  }
  ## Get starting values
  myparnames <- coef(x,mean=TRUE)
  paragree <- FALSE
  paragree.2 <- c()
  if (!is.null(optim$start)) {
    paragree <- myparnames%in%names(optim$start)
    paragree.2 <- names(optim$start)%in%myparnames
  }
  if (sum(paragree)>=length(myparnames))
    optim$start <- optim$start[which(paragree.2)]

 
  if (! (length(optim$start)==length(myparnames) & sum(paragree)==0)) 
  if (is.null(optim$start) || sum(paragree)<length(myparnames)) {
    Debug(list("starter=",optim$starterfun))
    start <- suppressWarnings(do.call(optim$starterfun, list(x=x,S=S,mu=mu,debug=lava.options()$debug,silent=silent)))
    Debug(start)  
    Debug(list("start=",start))
    if (length(paragree.2)>0) {
      start[which(paragree)] <- optim$start[which(paragree.2)]
    }
    optim$start <- start
  }

  ## Missing data
  if (missing) {
    control$start <- optim$start
    return(estimate.MAR(x=x,data=data,fix=fix,control=control,debug=lava.options()$debug,silent=silent,estimator=estimator,weight=weight,weight2=weight2,cluster=cluster,...))
  }
  

  ## Setup optimization constraints
  lowmin <- -Inf
  lower <- rep(lowmin,length(optim$start))
  if (length(optim$constrain)==1 & optim$constrain)
    lower[variances(x)+index(x)$npar.mean] <- optim$tol
  if (any(optim$constrain)) {
    if (length(optim$constrain)!=length(lower))
      constrained <- is.finite(lower)
    else
      constrained <- optim$constrain
    lower[] <- -Inf
      optim$constrain <- TRUE
    nn <- names(optim$start)
    CS <- optim$start[constrained]
    CS[CS<0] <- 0.01
    optim$start[constrained] <- log(CS)
    names(optim$start) <- nn
  }
  ## Fix problems with starting values? 
  optim$start[is.nan(optim$start)] <- 0  
  Debug(list("lower=",lower))
  
  ObjectiveFun  <- paste(estimator, "_objective", ".lvm", sep="")
  GradFun  <- paste(estimator, "_gradient", ".lvm", sep="")
  if (!exists(ObjectiveFun) & !exists(GradFun)) stop("Unknown estimator.")

  InformationFun <- paste(estimator, "_hessian", ".lvm", sep="")

  mymodel <- x  
  myclass <- "lvmfit"

  ## Non-linear parameter constraints involving observed variables? (e.g. nonlinear regression)
  xconstrain <- intersect(unlist(lapply(constrain(x),function(z) attributes(z)$args)),manifest(x))

  ## Random slopes or non-linear constraints?
  if (length(xfix)>0 | length(xconstrain)>0) { ## Yes
    x0 <- x
    
    if (length(xfix)>0) {
      myclass <- c("lvmfit.randomslope",myclass)
      nrow <- length(vars(x))
      xpos <- lapply(xfix,function(y) which(regfix(x)$labels==y))
      colpos <- lapply(xpos, function(y) ceiling(y/nrow))
      rowpos <- lapply(xpos, function(y) (y-1)%%nrow+1)
      myfix <- list(var=xfix, col=colpos, row=rowpos)
      x0 <- x
      for (i in 1:length(myfix$var)) 
        for (j in 1:length(myfix$col[[i]])) 
          regfix(x0, from=vars(x0)[myfix$row[[i]][j]],
                 to=vars(x0)[myfix$col[[i]][j]]) <-
                   colMeans(data[,myfix$var[[i]],drop=FALSE])
      x0 <- updatelvm(x0,zeroones=TRUE,deriv=TRUE)
      x <- x0
      yvars <- endogenous(x0)
      
      ## Alter start-values/constraints:
      new.par.idx <- coef(mymodel,mean=TRUE)%in%coef(x0,mean=TRUE)
      if (length(optim$start)>sum(new.par.idx))
        optim$start <- optim$start[new.par.idx]
      lower <- lower[new.par.idx]
      if (optim$constrain) {
        constrained <- constrained[new.par.idx]
      }
    }
    
    mydata <- as.matrix(data[,manifest(x0)])

    myObj <- function(pp) {
      if (optim$constrain) {
        pp[constrained] <- exp(pp[constrained])
      }
      myfun <- function(ii) {
        if (length(xfix)>0)
        for (i in 1:length(myfix$var)) {          
            x0$fix[cbind(rowpos[[i]],colpos[[i]])] <- index(x0)$A[cbind(rowpos[[i]],colpos[[i]])] <- data[ii,xfix[i]]
        }
        if (is.list(weight2)) {
          res <- do.call(ObjectiveFun, list(x=x0, p=pp, data=mydata[ii,], n=1, weight=weight[ii,], weight2=weight2[ii,]))
        } else
        {
          res <- do.call(ObjectiveFun, list(x=x0, p=pp, data=mydata[ii,], n=1, weight=weight[ii,], weight2=weight2))
        }
        return(res)
      }           
      sum(sapply(1:nrow(data),myfun))
    }

    myGrad <- function(pp) {
      if (optim$constrain) {
        pp[constrained] <- exp(pp[constrained])
      }
      myfun <- function(ii) {
        if (length(xfix)>0)
        for (i in 1:length(myfix$var)) {
          x0$fix[cbind(rowpos[[i]],colpos[[i]])] <- index(x0)$A[cbind(rowpos[[i]],colpos[[i]])] <- data[ii,xfix[i]]
        }
        if (is.list(weight2)) {
          rr <- do.call(GradFun, list(x=x0, p=pp, data=mydata[ii,,drop=FALSE], n=1, weight=weight[ii,], weight2=weight2))
        } else
        {
          rr <- do.call(GradFun, list(x=x0, p=pp, data=mydata[ii,,drop=FALSE], n=1, weight=weight[ii,], weight2=weight2[ii,]))
        }
        return(rr)        
      }
      ss <- rowSums(rbind(sapply(1:nrow(data),myfun)))
      if (optim$constrain) {
        ss[constrained] <- ss[constrained]*pp[constrained]
      }
      return(ss)
    }

    
    myInfo <- function(pp,...) {
      myfun <- function(ii) {
        if (length(xfix)>0)
        for (i in 1:length(myfix$var)) {
          x0$fix[cbind(rowpos[[i]],colpos[[i]])] <- index(x0)$A[cbind(rowpos[[i]],colpos[[i]])] <- data[ii,xfix[i]]
        }
        if (is.list(weight2)) {
          res <- do.call(InformationFun, list(p=pp, obj=myObj, x=x0, data=data[ii,],
                                              n=1, weight=weight[ii,], weight2=weight2))
        } else {
          res <- do.call(InformationFun, list(p=pp, obj=myObj, x=x0, data=data[ii,],
                                              n=1, weight=weight[ii,], weight2=weight2[ii,]))
        }
        return(res)
      }      
      L <- lapply(1:nrow(data),function(x) myfun(x))
      val <- apply(array(unlist(L),dim=c(length(pp),length(pp),nrow(data))),c(1,2),sum)
      if (!is.null(attributes(L[[1]])$grad)) {
        attributes(val)$grad <- colSums (
                                          matrix( unlist(lapply(L,function(i) attributes(i)$grad)) , ncol=length(pp), byrow=TRUE)
                                          )        
      }     
      return(val)
    }
    
##################################################
  } else { ## No, standard model

    myObj <- function(pp) {
      if (optim$constrain) {
        pp[constrained] <- exp(pp[constrained])
      }
      do.call(ObjectiveFun, list(x=x, p=pp, data=data, S=S, mu=mu, n=n, weight=weight
                                 ,weight2=weight2
                                 ))
    }
    myGrad <- function(pp) {
      if (optim$constrain)
        pp[constrained] <- exp(pp[constrained])
      S <- do.call(GradFun, list(x=x, p=pp, data=data, S=S, mu=mu, n=n, weight=weight
                                 , weight2=weight2
                                 ))
      if (optim$constrain) {
        S[constrained] <- S[constrained]*pp[constrained]
      }
      if (is.null(mu) & index(x)$npar.mean>0) {
        return(S[-c(1:index(x)$npar.mean)])
      }      
      return(S)
    }
    myInfo <- function(pp,...) {
      I <- do.call(InformationFun, list(p=pp, obj=myObj, x=x, data=data,
                                        S=S, mu=mu, n=n,
                                        weight=weight, weight2=weight2,
                                        type=optim$information))
      if (is.null(mu) & index(x)$npar.mean>0) {
        return(I[-c(1:index(x)$npar.mean),-c(1:index(x)$npar.mean)])
      }
      return(I)
    }
    
  }

  ## if (!exists(GradFun) & !is.null(optim$method)) {
  ##   message("Using numerical derivatives...\n")
  ##   myGrad <- function(pp) {
  ##     if (optim$constrain)
  ##       pp[constrained] <- exp(pp[constrained])
  ##     if (!require("numDeriv")) {        
  ##       S <- naiveGrad(myObj, pp)
  ##     } else {
  ##       S <- grad(myObj, pp, method=lava.options()$Dmethod)
  ##     }
  ##     if (optim$constrain) {
  ##       S[constrained] <- S[constrained]*pp[constrained]
  ##     }
  ##     return(S)        
  ##   }
  ##  }
##   if (!exists(InformationFun) & !is.null(optim$method)) {
##     if (!require("numDeriv")) stop("I do not know how to calculate the asymptotic variance of this estimator.
## For numerical approximation please install the library 'numDeriv'.")
##     message("Using a numerical approximation of hessian...\n");
##     myInfo <- function(pp,...) hessian(myObj, opt$estimate, method="Richardson")
##   }
  
  myHess <- function(pp) {
    p0 <- pp
    if (optim$constrain)
      pp[constrained] <- exp(pp[constrained])
    I0 <- myInfo(pp)
    attributes(I0)$grad <- NULL
    D <- attributes(I0)$grad
    if (is.null(D)) {
      D <- myGrad(p0)
      attributes(I0)$grad <- D
    }
    if (optim$constrain) {
      I0[constrained,-constrained] <- apply(I0[constrained,-constrained,drop=FALSE],2,function(x) x*pp[constrained]);
      I0[-constrained,constrained] <- t(I0[constrained,-constrained])
      if (sum(constrained)==1) {
        I0[constrained,constrained] <- I0[constrained,constrained]*outer(pp[constrained],pp[constrained]) - D[constrained]
      } else {
        I0[constrained,constrained] <- I0[constrained,constrained]*outer(pp[constrained],pp[constrained]) - diag(D[constrained])
      }
    }
    return(I0)
  }   
  
##  if (!exists(InformationFun)) myInfo <- myHess <- NULL
##  if (is.null(get(InformationFun))) myInfo <- myHess <- NULL
  if (is.null(tryCatch(get(InformationFun),error = function (x) NULL)))
    myInfo <- myHess <- NULL
  if (is.null(tryCatch(get(GradFun),error = function (x) NULL)))
    myGrad <- NULL

  coefname <- coef(x,mean=optim$meanstructure);
  ##  browser()
  if (!silent) message("Optimizing objective function...")
  if (optim$trace>0 & !silent) message("\n")
  ## Optimize with lower constraints on the variance-parameters
  if ((is.data.frame(data) | is.matrix(data)) && nrow(data)==0) stop("No observations")

  if (!is.null(optim$method)) {
    opt <- do.call(optim$method,
                   list(start=optim$start, objective=myObj, gradient=myGrad, hessian=myHess, lower=lower, control=optim, debug=debug))
    if (is.null(opt$estimate))
      opt$estimate <- opt$par
    if (optim$constrain) {
      opt$estimate[constrained] <- exp(opt$estimate[constrained])
    }
    names(opt$estimate) <- coefname
    
    opt$gradient <- as.vector(myGrad(opt$par))
  } else {
    opt <- do.call(ObjectiveFun, list(x=x,data=data,control=control,...))
    opt$grad <- rep(0,length(opt$estimate))
  }
  if (quick) {
    return(opt$estimate)
  }
  ## Calculate std.err:
  
  pp <- rep(NA,length(coefname)); names(pp) <- coefname
  pp[names(opt$estimate)] <- opt$estimate
  pp.idx <- na.omit(match(coefname,names(opt$estimate)))

  mom <- tryCatch(modelVar(x, pp, data=data),error=function(x)NULL)
  if (!silent) message("\nCalculating asymptotic variance...\n")
  asVarFun  <- paste(estimator, "_variance", ".lvm", sep="")
  if (!exists(asVarFun)) {
    if (is.null(myInfo)) {
      if (!is.null(myGrad))
        myInfo <- function(pp,...)
          jacobian(myGrad,pp,method=lava.options()$Dmethod)
      else
        myInfo <- function(pp,...)
          -hessian(myObj,pp,method=lava.options()$Dmethod)
    }
    I <- myInfo(opt$estimate)
    asVar <- tryCatch(solve(I),
                      error=function(e) matrix(NA, length(opt$estimate), length(opt$estimate)))
##    diag(asVar)[(diag(asVar)==0)] <- NA
  } else {
    asVar <- tryCatch(do.call(asVarFun,
                              list(x=x,p=opt$estimate,data=data,opt=opt)),
                      error=function(e) matrix(NA, length(opt$estimate), length(opt$estimate)))
    
  }

  if (any(is.na(asVar))) {warning("Problems with asymptotic variance matrix. Possibly non-singular information matrix!")
                        }
  diag(asVar)[(diag(asVar)==0)] <- NA

  Debug("did that") 

  nparall <- index(x)$npar + ifelse(optim$meanstructure, index(x)$npar.mean,0)
  mycoef <- matrix(NA,nrow=nparall,ncol=4)

  mycoef[pp.idx,1] <- opt$estimate
  
  ### OBS: v = t(A)%*%v + e
  res <- list(model=x, call=cl, coef=mycoef, vcov=asVar, mu=mu, S=S, ##A=A, P=P,
              model0=mymodel, ## Random slope hack
              estimator=estimator, opt=opt,
              data=list(model.frame=data, S=S, mu=mu, C=mom$C, v=mom$v, n=n,
                m=length(latent(x)), k=k),
              weight=weight, weight2=weight2,
              cluster=cluster,
              pp.idx=pp.idx,
              graph=NULL, control=optim)

  class(res) <- myclass
  
  myhooks <- gethook("post.hooks")
  for (f in myhooks) {
    res0 <- do.call(f,list(x=res))
    if (!is.null(res0))
      res <- res0
  }
  
  if(graph) {
    res <- edgelabels(res,type="est")
  }
  
  return(res)
}

###}}} estimate.lvm

###{{{ estimate.formula

##' @S3method estimate formula
estimate.formula <- function(x,data,pred.norm=c(),unstruct=FALSE,silent=TRUE,...) {
  cl <- match.call()
  ## {  varnames <- all.vars(x)
  ##    mf <- model.frame(x,data)
  ##    mt <- attr(mf, "terms")
  ##    yvar <- names(mf)[1]
  ##    y <- data[,yvar]
  ##    opt <- options(na.action="na.pass")
  ##    mm <- model.matrix(x,data)
  ##    options(opt)
  ##    covars <- colnames(mm)
  ##    if (attr(terms(x),"intercept")==1)
  ##      covars <- covars[-1]
  ##    model <- lvm()
  ##    for (i in covars) {
  ##      model <- regression(model, to=yvar, from=i,silent=TRUE)
  ##    }     
  ##    mydata <- as.data.frame(cbind(y,mm)); names(mydata)[1] <- yvar
  ##  }
  model <- lvm(x,silent=silent)
  ##  covars <- exogenous(model)
  ##  exogenous(model) <- setdiff(covars,pred.norm)
  ##  if (unstruct) {    
  ##    model <- covariance(model,pred.norm,pairwise=TRUE)
  ##  }
  estimate(model,data,silent=silent,...)
}

###}}} estimate.formula

###{{{ estimate.lm

##' @S3method estimate lm
estimate.lm <- function(x,data=model.frame(x),...) {
  estimate(formula(x),data,...)
}

###}}} estimate.lm

