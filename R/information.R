`information` <-
function(x,...) UseMethod("information")

###{{{ information.lvm

information.lvm <- function(x,p,n,type=ifelse(model=="gaussian",
                                    c("E","hessian","varS","outer","robust","num"),"outer"),
                            data,weight=NULL,model="gaussian",method="simple",
                            score=TRUE,...) {
  if (missing(n))
    n <- NROW(data)
  if (type[1]=="robust") {
    cl <- match.call()
    cl$type <- "E"
    I <- eval.parent(cl)
    cl$type <- "outer"
    J <- eval.parent(cl)
##    I <- information(x,p,n,data=data,weight=weight,type="E")
##    J <- information(x,p,n,data=data,weight=weight,type="outer")
    return(I%*%solve(J)%*%I)
  }
  if (type[1]%in%c("num","hessian")  | (type[1]%in%c("E","hessian") & model!="gaussian")) {
    require("numDeriv")
    myf <- function(p0) score(x, p=p0, model=model,data=data, weight=weight,indiv=FALSE,seed=1,n=n) ##...)
    ##    I <- -hessian(function(p0) logLik(x,p0,dd),p)
    I <- -jacobian(myf,p,method=method)
    return((I+t(I))/2) # Symmetric result
  }
  if (type[1]=="varS" | type[1]=="outer") {
    S <- score(x,p=p,data=na.omit(data),model=model,weight=weight,indiv=TRUE,...)
##    print("...")
    res <- t(S)%*%S
    attributes(res)$grad <- colSums(S)
    return(res)    
  }


  if (n>1) {
    xfix <- colnames(data)[(colnames(data)%in%parlabels(x))]
    xconstrain <- intersect(unlist(lapply(constrain(x),function(z) attributes(z)$args)),manifest(x))

    if (length(xfix)>0 | length(xconstrain)>0) { ##### Random slopes!
      x0 <- x
      if (length(xfix)>0) {
        nrow <- length(vars(x))
        xpos <- lapply(xfix,function(y) which(regfix(x)$labels==y))
        colpos <- lapply(xpos, function(y) ceiling(y/nrow))
        rowpos <- lapply(xpos, function(y) (y-1)%%nrow+1)
        myfix <- list(var=xfix, col=colpos, row=rowpos)
        for (i in 1:length(myfix$var)) 
          for (j in 1:length(myfix$col[[i]])) 
            regfix(x0, from=vars(x0)[myfix$row[[i]]][j],to=vars(x0)[myfix$col[[i]]][j]) <-
              data[1,myfix$var[[i]]]
        ##rep(data[1,myfix$var[[i]]],length(myfix$row[[i]]))
        index(x0) <- reindex(x0,zeroones=TRUE,deriv=TRUE)
      }
      pp <- modelPar(x0,p)
      p0 <- with(pp, c(meanpar,p))
      k <- length(index(x)$manifest)      
      myfun <- function(ii) {
        if (length(xfix)>0)
          for (i in 1:length(myfix$var)) {
            for (j in 1:length(myfix$col[[i]])) {
              index(x0)$A[cbind(myfix$row[[i]],myfix$col[[i]])] <- data[ii,myfix$var[[i]]]
            }
          }
        ww <- NULL
        if (!is.null(weight))
          ww <- weight[ii,]
        return(information(x0,p=p,n=1,type=type,weight=ww,data=data[ii,]))
      }
      L <- lapply(1:nrow(data),function(y) myfun(y))
      val <- apply(array(unlist(L),dim=c(length(p0),length(p0),nrow(data))),c(1,2),sum)
      return(val)
    }    
  }

  if (!is.null(weight) && is.matrix(weight)) {
    L <- lapply(1:nrow(weight),function(y) information(x,p=p,n=1,type=type,weight=weight[y,]))
    val <- apply(array(unlist(L),dim=c(length(p),length(p),nrow(weight))),c(1,2),sum)
    return(val)
  }
  
  mp <- moments(x,p,data=data)
  pp <- modelPar(x,p)
  D <- deriv(x, meanpar=pp$meanpar, mom=mp, p=p)##, all=length(constrain(x))>0)
  C <- mp$C
  iC <- Inverse(C,0,det=FALSE)

  if (is.null(weight)) {
    W <- diag(ncol(iC))
  } else {
    if (length(weight)<ncol(iC)) {
      oldweight <- weight
      weight <- rbind(rep(1,ncol(iC))) ## Ones at exogenous var.
      idx <- index(x)$vars%in%index(x)$exogenous
      print(idx); print(oldweight)
      weight[,idx] <- oldweight
    }
    W <- diag(as.numeric(weight))
    iW <- W
    diag(iW) <- 1/diag(iW)
  }


  if (NCOL(D$dS)>30) {
    I <- matrix(0,NCOL(D$dS),NCOL(D$dS))
    for (i in 1:NCOL(D$dS)) {
      for (j in i:NCOL(D$dS)) {        
        I[i,j] <- I[j,i] <- 
          sum(diag(matrix(D$dS[,i],NCOL(iC))%*%iC%*%W%*%matrix(D$dS[,j],NCOL(iC))%*%iC))
      }
    }
    information_Sigma <- n/2*I
      
  } else {
    information_Sigma <-  n/2*t(D$dS)%*%((iC)%x%(iC%*%W))%*%(D$dS)
  }
  if (is.null(pp$meanpar))
    return(information_Sigma)

  f <- function(p0) modelVar(x,p0)$xi
  
  idx <- with(index(x), (1:npar) + npar.mean)
  dxi <- D$dxi; ##dxi[,idx] <- 0
  information_mu <- n*t(D$dxi) %*% (iC%*%W) %*% (D$dxi)
  
  information <- information_Sigma + information_mu
  ##browser()
  ## if (score) {
  ##   vec.iC <- as.vector(iC)  
  ##   Grad <- n/2*crossprod(D$dS, as.vector(iC%*%T%*%iC)-vec.iC)
  ##   if (!is.null(mu)) # & mp$npar.mean>0)
  ##     Grad <- Grad - n/2*crossprod(D$dT,vec.iC)
  ## }

  return(information)  
}

###}}} information.lvm

###{{{ information.lvmfit

information.lvmfit <- function(x,p=pars(x),n=x$data$n,data=model.frame(x),model=x$estimator,weight=Weight(x),...) {
  information(x$model0,p=p,n=n,data=data,model=model,weight=weight,...)
}

###}}} information.lvmfit

