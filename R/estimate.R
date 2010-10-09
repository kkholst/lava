`estimate` <-
estimate <- function(x,...) UseMethod("estimate")

###{{{ estimate.lvm

`estimate.lvm` <-
function(x, data,
         estimator="gaussian",
         control=list(),
         missing=FALSE,
         weight,
         index=TRUE,
         graph=FALSE,
         fix,
         debug=FALSE, silent=FALSE,
         quick=FALSE,
         ...) { 

  if (length(exogenous(x)>0)) {
    catx <- categorical2dummy(x,data)
    x <- catx$x; data <- catx$data
  }
  
  cl <- match.call()
  Method <-  paste(estimator, "_method", ".lvm", sep="")
  if (!exists(Method))
    Method <- "nlminb2"
  else
    Method <- get(Method)
  optim <- list(
                eval.max=300,
                iter.max=500,
                trace=ifelse(debug,3,0),
                gamma=1,
                gamma2=1,
                ngamma=NULL,
                lambda=0.05,
                abs.tol=1e-12,
                epsilon=1e-10,
                delta=1e-10,
                S.tol=1e-5,
                stabil=FALSE,
                start=NULL,
                constrain=FALSE,
                method=Method,
                starterfun="startvalues",
                information="E",
                meanstructure=TRUE,
                sparse=FALSE,
                tol=1e-9)

  if (length(control)>0) {
    optim[names(control)] <- control
  }
  ##  if (length(optcontrol)>0) {
  ##    optim$control[names(optcontrol)] <- optcontrol
  ##  }  
  ## Random-slopes:
  redvar <- intersect(intersect(parlabels(x),latent(x)),colnames(data))
  if (length(redvar)>0)
    warning(paste("Remove latent variable colnames from dataset",redvar))
  xfix <- setdiff(colnames(data)[(colnames(data)%in%parlabels(x,exo=TRUE))],latent(x))
  if (missing(fix)) {
    fix <- ifelse(length(xfix)>0,FALSE,TRUE)
  }
  Debug(list("start=",optim$start),debug)
  ## Weights...
  if (!missing(weight)) {
    if (is.character(weight)) {
      weight <- data[,weight]
    }
    weight <- cbind(weight)
  } else {
    weight <- NULL
  }
  Debug("procdata",debug)
  dd <- procdata.lvm(x,data=data,debug=debug)
  S <- dd$S; mu <- dd$mu; n <- dd$n
  Debug(list("n=",n),debug)  
  Debug(list("S=",S),debug)
  Debug(list("mu=",mu),debug)  
  ### Run hooks (additional lava plugins)
  myhooks <- gethook()
  for (f in myhooks) {
    res <- do.call(f, list(x=x,data=data,weight=weight,estimator=estimator,optim=optim))
    if (!is.null(res$x)) x <- res$x
    if (!is.null(res$data)) data <- res$data
    if (!is.null(res$weight)) weight <- res$weight
    if (!is.null(res$optim)) optim <- res$optim
    if (!is.null(res$estimator)) estimator <- res$estimator
  }
  
  if (!quick & index) {
    ## Proces data and setup some matrices
    x <- fixsome(x, measurement.fix=fix, S=S, mu=mu, n=n,debug=!silent)
    if (!silent)
    cat("Reindexing model...\n")
    if (length(xfix)>0) {
      index(x) <- reindex(x,sparse=optim$sparse,zeroones=TRUE,deriv=TRUE,debug=debug)
    } else {
      x <- updatelvm(x,sparse=optim$sparse,zeroones=TRUE,deriv=TRUE,mean=TRUE,debug=debug)
    }
  }
  
##  if (index(x)$npar.mean==0) optim$meanstructure <- FALSE       
  if (is.null(estimator) || estimator==FALSE) {
    return(x)
  }  
  k <- length(manifest(x))
  ##  m <- length(latent(x))    
  Debug(list("S=",S),debug)
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
    Debug(list("starter=",optim$starterfun), debug)
    start <- suppressWarnings(do.call(optim$starterfun, list(x=x,S=S,mu=mu,debug=debug,silent=silent)))
    Debug(start, debug)  
    Debug(list("start=",start),debug)
    if (length(paragree.2)>0) {
      start[which(paragree)] <- optim$start[which(paragree.2)]
    }
    optim$start <- start
  }
  
  ## Missing data
  if (missing) {
    control$start <- optim$start
    return(estimate.MAR(x=x,data=data,fix=fix,control=control,debug=debug,silent=silent,estimator=estimator,weight=weight,...))
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
    CS <- optim$start[constrained]
    CS[CS<0] <- 0.01
    optim$start[constrained] <- log(CS)
  }
  ## Fix problems with starting values? 
  optim$start[is.nan(optim$start)] <- 0  
  Debug(list("lower=",lower),debug)

  
  ObjectiveFun  <- paste(estimator, "_objective", ".lvm", sep="")
  if (!exists(ObjectiveFun)) stop("Unknown estimator.")
  GradFun  <- paste(estimator, "_gradient", ".lvm", sep="")
  InformationFun <- paste(estimator, "_hessian", ".lvm", sep="")

  mymodel <- x  
  myclass <- "lvmfit"

  ## Non-linear parameter constraints involving observed variables? (e.g. nonlinear regression)
  xconstrain <- intersect(unlist(lapply(constrain(x),function(z) attributes(z)$args)),manifest(x))
  
###### Random slopes or non-linear constraints?
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
        res <- do.call(ObjectiveFun, list(x=x0, p=pp, data=mydata[ii,], n=1, weight=weight[ii,]))
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
        rr <- do.call(GradFun, list(x=x0, p=pp, data=mydata[ii,,drop=FALSE], n=1, weight=weight[ii,]))
##        browser()

        ##       rr <- score(x0,p=pp,data=data[ii,])
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
##        browser()
        res <- do.call(InformationFun, list(p=pp, obj=myObj, model=x0, data=data[ii,],
                                            n=1, weight=weight[ii,]))
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
      do.call(ObjectiveFun, list(x=x, p=pp, data=data, S=S, mu=mu, n=n, weight=weight))
    }
    myGrad <- function(pp) {
      if (optim$constrain)
        pp[constrained] <- exp(pp[constrained])
      S <- do.call(GradFun, list(x=x, p=pp, data=data, S=S, mu=mu, n=n, weight=weight))
      if (optim$constrain) {
        S[constrained] <- S[constrained]*pp[constrained]
      }
      if (is.null(mu) & index(x)$npar.mean>0) {
        return(S[-c(1:index(x)$npar.mean)])
      }      
      return(S)
    }
    myInfo <- function(pp,...) {
      I <- do.call(InformationFun, list(p=pp, obj=myObj, model=x, data=data,
                                        S=S, mu=mu, n=n,
                                        weight=weight, type=optim$information))
      if (is.null(mu) & index(x)$npar.mean>0) {
        return(I[-c(1:index(x)$npar.mean),-c(1:index(x)$npar.mean)])
      }
      return(I)
    }
    
  }
  if (!exists(GradFun) & !is.null(optim$method)) {
    cat("Using numerical derivatives...\n")
    myGrad <- function(pp) {
      if (optim$constrain)
        pp[constrained] <- exp(pp[constrained])
      if (!require("numDeriv")) {        
        S <- naiveGrad(myObj, pp)
      } else {
        S <- grad(myObj, pp, method="Richardson")
      }
      if (optim$constrain) {
        S[constrained] <- S[constrained]*pp[constrained]
      }
      return(S)        
    }
  }  
  if (!exists(InformationFun) & !is.null(optim$method)) {
    if (!require("numDeriv")) stop("I do not know how to calculate the asymptotic variance of this estimator.
For numerical approximation please install the library 'numDeriv'.")
    cat("Using a numerical approximation of hessian...\n");
    myInfo <- function(pp,...) hessian(myObj, opt$estimate, method="Richardson")
  }
  
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
  coefname <- coef(x,mean=optim$meanstructure);

  if (!silent) cat("Optimizing objective function...\n")
  ## Optimize with lower constraints on the variance-parameters
  if (!is.null(optim$method)) {
    opt <- do.call(optim$method,
                   list(start=optim$start, objective=myObj, gradient=myGrad, hessian=myHess, lower=lower, control=optim, debug=debug))
    opt$estimate <- opt$par
    if (optim$constrain) {
      opt$estimate[constrained] <- exp(opt$estimate[constrained])
    }
    names(opt$estimate) <- coefname 
    opt$gradient <- as.vector(myGrad(opt$par))
  } else {
    opt <- do.call(ObjectiveFun, list(x=x,data=data))
    opt$grad <- rep(0,length(opt$estimate))
  }
  if (quick) return(opt$estimate)
  ## Calculate std.err:
  
  pp <- rep(NA,length(coefname)); names(pp) <- coefname
  pp[names(opt$estimate)] <- opt$estimate
  pp.idx <- na.omit(match(coefname,names(opt$estimate)))

  mom <- modelVar(x, pp, debug)
  if (!silent) cat("Calculating asymptotic variance...\n")
  asVarFun  <- paste(estimator, "_variance", ".lvm", sep="")
  if (!exists(asVarFun)) {
    asVar <- Inverse(myInfo(opt$estimate))
    diag(asVar)[(diag(asVar)==0)] <- NA
##    if ()
##    asVar <- tryCatch(
##                      solve(myInfo(opt$estimate)),
##                      error=function(e) matrix(NA, length(opt$estimate), length(opt$estimate)))
###
  } else {
    asVar <- tryCatch(do.call(asVarFun,
                              list(x=x,p=opt$estimate,data=data,opt=opt)),
                      error=function(e) matrix(NA, length(opt$estimate), length(opt$estimate)))
    
  }
  if (any(is.na(asVar))) {warning("Problems with asymptotic variance matrix. Possibly non-singular information matrix!")
                        }
  Debug("did that", debug) 

  ##
  SD <- sqrt(diag(asVar))
  Z <- opt$estimate/SD
  pval <- 2*(1-pnorm(abs(Z)))
  coef <- cbind(opt$estimate, SD, Z, pval);
  Debug(coef,debug)

  nparall <- mom$npar + ifelse(optim$meanstructure, mom$npar.mean,0)
  mycoef <- matrix(NA,nrow=nparall,ncol=4)
  mycoef[pp.idx,] <- coef
  colnames(mycoef) <- c("Estimate","Std. Error", "Z value", "Pr(>|z|)")    

  mynames <- c()
  if (optim$meanstructure) {
    mynames <- vars(x)[index(x)$v1==1]
  }
  if (index(x)$npar>0) {
    mynames <- c(mynames,paste("p",1:index(x)$npar,sep=""))
  }
  rownames(mycoef) <- mynames
  Debug(list("COEF",coef),debug)
   
  ### OBS: v = t(A)%*%v + e
  res <- list(model=x, call=cl, coef=mycoef, vcov=asVar, mu=mu, S=S, ##A=A, P=P,
              model0=mymodel, ## Random slope hack
              estimator=estimator, opt=opt,
              data=list(model.frame=data, S=S, mu=mu, C=mom$C, v=mom$v, n=n,
                m=length(latent(x)), k=k),
              weight=weight,
              graph=NULL, control=optim)

  class(res) <- myclass
  if(graph) {
    res <- edgelabels(res,type="est")
  }
  return(res)
}

###}}} estimate.lvm

###{{{ estimate.formula

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
  model <- lvm()
  regression(model,silent=silent) <- x
  ##  covars <- exogenous(model)
  ##  exogenous(model) <- setdiff(covars,pred.norm)
  ##  if (unstruct) {    
  ##    model <- covariance(model,pred.norm,pairwise=TRUE)
  ##  }
  estimate(model,data,silent=silent,...)
}

###}}} estimate.formula

###{{{ estimate.lm

estimate.lm <- function(x,data=model.frame(x),...) {
  estimate(formula(x),data,...)
}

###}}} estimate.lm

