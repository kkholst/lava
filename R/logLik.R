###{{{ logLik.lvm

logLik.lvm <- function(object,p,data,model="gaussian",indiv=FALSE,S,mu,n,debug=FALSE,weight=NULL,...) {
  cl <- match.call()
  xfix <- colnames(data)[(colnames(data)%in%parlabels(object,exo=TRUE))]
  xconstrain <- intersect(unlist(lapply(constrain(object),function(z) attributes(z)$args)),manifest(object))


  
  Debug(xfix,debug)
  if (missing(n)) {
    n <- nrow(data)
  }
  lname <- paste(model,"_logLik.lvm",sep="")
  logLikFun <- get(lname)
  if (length(xfix)>0 | length(xconstrain)>0) { ##### Random slopes!
    x0 <- object
    if (length(xfix)>0) {
      Debug("random slopes...",debug)
      nrow <- length(vars(object))
      xpos <- lapply(xfix,function(y) which(regfix(object)$labels==y))
      colpos <- lapply(xpos, function(y) ceiling(y/nrow))
      rowpos <- lapply(xpos, function(y) (y-1)%%nrow+1)
      myfix <- list(var=xfix, col=colpos, row=rowpos)
      for (i in 1:length(myfix$var))
        ##      regfix(x0, from=vars(x0)[myfix$row[[i]][]],to=vars(x0)[myfix$col[[i]][j]]) <
        ##          (data[1,myfix$var[[i]]])
        for (j in 1:length(myfix$col[[i]])) {
          regfix(x0, from=vars(x0)[myfix$row[[i]][j]],to=vars(x0)[myfix$col[[i]][j]]) <-
            data[1,myfix$var[[i]]]
        }
      index(x0) <- reindex(x0,zeroones=TRUE,deriv=TRUE)
    }
    pp <- modelPar(x0,p)
    p0 <- with(pp, c(meanpar,p))
    k <- length(index(x0)$manifest)
    myfun <- function(ii) {
      if (length(xfix)>0)
        for (i in 1:length(myfix$var)) {
          index(x0)$A[cbind(myfix$row[[i]],myfix$col[[i]])] <- data[ii,myfix$var[[i]]]
        }
      return(logLikFun(x0,data=data[ii,,drop=FALSE], p=with(pp,c(meanpar,p)),weight=weight[ii,,drop=FALSE],model=model,debug=debug,indiv=indiv,...))
    }    
    loglik <- sapply(1:nrow(data),myfun)
    if (!indiv) {
      loglik <- sum(loglik)
      n <- nrow(data)
      attr(loglik, "nall") <- n
      attr(loglik, "nobs") <- n-length(p0)
      attr(loglik, "df") <- length(p0)
      class(loglik) <- "logLik"      
    }
    return(loglik)
  }

  cl[[1]] <- logLikFun
  loglik <- eval.parent(cl)

  if (is.null(attr(loglik,"nall")))
    attr(loglik, "nall") <- n
  if (is.null(attr(loglik,"nobs")))
    attr(loglik, "nobs") <- n-length(p)
  if (is.null(attr(loglik,"df")))
    attr(loglik, "df") <- length(p)
  class(loglik) <- "logLik"
  return(loglik)
}

###}}}

###{{{ gaussian_loglik

gaussian_logLik.lvm <- function(object,p,data,
                          type=c("cond","sim","exo","sat","cond2"),
                          weight=NULL, indiv=FALSE, S, mu, n, debug=FALSE, meanstructure=TRUE,...) { 
   
  exo.idx <- with(index(object), exo.obsidx)##match(exogenous(object),manifest(object))
  endo.idx <- with(index(object), endo.obsidx)##match(endogenous(object),manifest(object))
  if (type[1]=="exo") {
    if (length(exo.idx)==0 || is.na(exo.idx))
      return(0)
  }
  
  cl <- match.call()  
  if (type[1]=="cond") {
    cl$type <- "sim"
    L0 <- eval.parent(cl)
    cl$type <- "exo"
    L1 <- eval.parent(cl)
    loglik <- L0-L1
    return(loglik)
  }
    
  if (missing(n)) {
    if (is.vector(data)) n <- 1
    else n <- nrow(data)
  }
  k <- length(index(object)$manifest)


  if (type[1]=="sat") {
    if (missing(S)) {
      d0 <- procdata.lvm(object,data=data)
      S <- d0$S; mu <- d0$mu; n <- d0$n
      
    }
    if (missing(p)) p <- rep(1,length(coef(object)))
    L1 <- logLik(object,p,data,type="exo",meanstructure=meanstructure)
    ##    Sigma <- (n-1)/n*S ## ML = 1/n * sum((xi-Ex)^2)
    Sigma <- S
    loglik <- -(n*k)/2*log(2*base::pi) -n/2*(log(det(Sigma)) + k) - L1
    P <- length(endo.idx)
    k <- length(exo.idx)
    npar <- P*(1+(P-1)/2)
    if (meanstructure) npar <- npar+ (P*k + P)
    attr(loglik, "nall") <- n
    attr(loglik, "nobs") <- n-npar
    attr(loglik, "df") <- npar    
    class(loglik) <- "logLik"    
    return(loglik)
  }
  myidx <- switch(type[1],
                  sim =  1:length(index(object)$manifest),
                  cond = { endo.idx },
                  cond2 = { endo.idx },
                  exo =  { exo.idx } )

  mom <- moments(object, p, conditional=(type[1]=="cond2"), data=data)  
  C <- mom$C
  xi <- mom$xi
  if (type[1]=="exo") {
    C <- C[exo.idx,exo.idx,drop=FALSE]
    xi <- xi[exo.idx,drop=FALSE]
  }  
  Debug(list("C=",C),debug)
  k <- nrow(C)

  iC <- Inverse(C,det=TRUE)
  detC <- attributes(iC)$det
#  iC <- solve(C)
#  detC <- 1
  
  ##  detC <- det(C)            
  ##  iC <- try(solve(C), silent=TRUE)
  ##   if (detC<0 | inherits(iC, "try-error")) {
  ##     if (indiv)
  ##       return(rep(-.Machine$double.xmax,nrow(data)))
  ##     else
  ##       return(-.Machine$double.xmax)
  ##   }


  if (!is.null(weight)) {
    weight <- cbind(weight)
    K <- length(exo.idx)+length(endo.idx)
    if (ncol(weight)!=1 & ncol(weight)!=K) {
      w.temp <- weight
      weight <- matrix(1,nrow=nrow(weight),ncol=K)
      weight[,endo.idx] <- w.temp
    }
    if (type=="exo")
      weight <- NULL
  }
  
  notdatalist <- (!is.list(data) | is.data.frame(data))
  if (missing(n))
    if (!missing(data)) n <- NROW(data)
  if (!missing(n))
  if (notdatalist & (n<2 | indiv | !is.null(weight))) {
    if (n==1)
      data <- rbind(data) 
    res <- numeric(n)
    data <- data[,index(object)$manifest,drop=FALSE]
    loglik <- 0; 
    for (i in 1:n) {
      ti <- cbind(as.numeric(data[i,myidx]))
      if (meanstructure) {
        ti <- ti-xi
      }
      if (!is.null(weight)) {
        W <- diag(weight[i,])
        val <- -k/2*log(2*base::pi) -1/2*log(detC) - 1/2*(t(ti)%*%W)%*%iC%*%(ti)
      } else { 
        val <- -k/2*log(2*base::pi) -1/2*log(detC) - 1/2*t(ti)%*%iC%*%(ti)
      }
      if (indiv)
        res[i] <- val
      loglik <- loglik + val
    }
    if (indiv)
      return(res)    
  } else {
   if (missing(S)) {
      d0 <- procdata.lvm(object,data=data)
      S <- d0$S; mu <- d0$mu; n <- d0$n
    }
    S <- S[myidx,myidx,drop=FALSE]
    mu <- mu[myidx,drop=FALSE]  
    T <- S
    if (meanstructure) {
      W <- tcrossprod(mu-xi)
      T <- S+W
    }
    loglik <- -(n*k)/2*log(2*base::pi) -n/2*(log(detC) + tr(T%*%iC))
  }
  return(loglik)
}

###}}}  

###{{{ logLik.lvmfit

logLik.lvmfit <- function(object,
                          p=coef(object),
                          data=model.frame(object),
                          model=object$estimator,
                          weight=Weight(object),
##                          meanstructure=object$control$meanstructure,
                          ...) {
  logLikFun <- paste(model,"_logLik.lvm",sep="")
  if (!exists(logLikFun)) {
    model <- "gaussian"
  }
  l <- logLik.lvm(object$model0,p,data,model=model,weight=weight,##meanstructure=meanstructure,...)
                  ...)
  return(l)
}

###}}} logLik.lvmfit

###{{{ logLik.lvm.missing

logLik.lvm.missing <- function(object,
                               p=pars(object), model=object$estimator,
                               weight=Weight(object$estimate),                
                               ...) {
  logLik(object$estimate$model0, p=p, model=model, weight=weight, ...)
}

###}}}

###{{{ logLik.multigroup

logLik.multigroup <- function(object,p,data=object$data,weight=NULL,type=c("cond","sim","exo","sat"),...) {
  res <- procrandomslope(object)
  pp <- with(res, modelPar(model,p)$p) 
    
  if (type[1]=="sat") {
    n <- 0
    df <- 0
    loglik <- 0
    for (i in 1:object$ngroup) {
      m <- Model(object)[[i]]
      L <- logLik(m,p=pp[[i]],data=object$data[[i]],type="sat")
      df <- df + attributes(L)$df
      loglik <- loglik + L
      n <- n + object$samplestat[[i]]$n
    }
    attr(loglik, "nall") <- n
    attr(loglik, "nobs") <- n-df
    attr(loglik, "df") <- df
    class(loglik) <- "logLik"
    return(loglik)  
  }

  n <- 0
  loglik <- 0; for (i in 1:object$ngroup) {
    n <- n + object$samplestat[[i]]$n
    val <- logLik(object$lvm[[i]],pp[[i]],data[[i]],weight=weight[[i]],type=type,...)
    loglik <- loglik + val
  }
  attr(loglik, "nall") <- n
  attr(loglik, "nobs") <- n-length(p)
  attr(loglik, "df") <- length(p)
  class(loglik) <- "logLik"
  return(loglik)  
}

###}}} logLik.multigroup

###{{{ logLik.multigroupfit
logLik.multigroupfit <- function(object,
                                 p=pars(object), weight=Weight(object), model=object$estimator, ...) {
  logLik(object$model0,p=p,weight=weight,model=model,...)
}
###}}} logLik.multigroup
