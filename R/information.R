`information` <-
function(x,...) UseMethod("information")

###{{{ information.lvm

information.lvm <- function(x,p,n,type=c("E","obs","varS","robust","num"),data,weight=NULL,...) {
  if (missing(n))
    n <- NROW(data)
  if (type[1]=="robust") {
    I <- information(x,p,n,data=data,weight=weight,type="E")
    J <- information(x,p,n,data=data,weight=weight,type="varS")
    return(I%*%solve(J)%*%I)
  }
  if (type[1]=="num") {
    require("numDeriv")
    myf <- function(p0) logLik(x, p0, ...)
    ##    I <- -hessian(function(p0) logLik(x,p0,dd),p)
    return(-hessian(myf,p))
  }
  if (type[1]=="varS") {
    S <- score(x,p=p,data=na.omit(data))
    return(t(S)%*%S)    
  }   
  
  if (!missing(data)) {
    xfix <- colnames(data)[(colnames(data)%in%parlabels(x))]
    if (length(xfix)>0) { ##### Random slopes!
      nrow <- length(vars(x))
      xpos <- lapply(xfix,function(y) which(regfix(x)$labels==y))
      colpos <- lapply(xpos, function(y) ceiling(y/nrow))
      rowpos <- lapply(xpos, function(y) (y-1)%%nrow+1)
      myfix <- list(var=xfix, col=colpos, row=rowpos)
      x0 <- x
      for (i in 1:length(myfix$var)) 
        for (j in 1:length(myfix$col[[i]])) 
          regfix(x0, from=vars(x0)[myfix$row[[i]]][j],to=vars(x0)[myfix$col[[i]]][j]) <-
            data[1,myfix$var[[i]]]
            ##rep(data[1,myfix$var[[i]]],length(myfix$row[[i]]))
      index(x0) <- reindex(x0,zeroones=TRUE,deriv=TRUE)
      pp <- modelPar(x0,p)
      p0 <- with(pp, c(meanpar,p))
      k <- length(index(x)$manifest)
      myfun <- function(ii) {
        for (i in 1:length(myfix$var)) {
          for (j in 1:length(myfix$col[[i]])) {
            index(x0)$A[cbind(myfix$row[[i]],myfix$col[[i]])] <- data[ii,myfix$var[[i]]]
          }
        }
        ww <- NULL
        if (!is.null(weight))
          ww <- weight[ii,]
        return(information(x0,p=p,n=1,type=type,weight=ww))
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
      
  information_Sigma <-  n/2*t(D$dS)%*%((iC)%x%(iC%*%W))%*%(D$dS)
  if (is.null(pp$meanpar))
    return(information_Sigma)

  f <- function(p0) modelVar(x,p0)$xi
  
  idx <- with(index(x), (1:npar) + npar.mean)
  dxi <- D$dxi; ##dxi[,idx] <- 0
  information_mu <- n*t(D$dxi) %*% (iC%*%W) %*% (D$dxi)

  information <- information_Sigma + information_mu
  return(information)  
}

###}}} information.lvm

###{{{ information.lvmfit

information.lvmfit <- function(x,p=pars(x),n=x$data$n,data=model.frame(x),...) {
  information(x$model0,p=p,n=n,data=data,...)
}

###}}} information.lvmfit

