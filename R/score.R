`score` <-
function(x,...) UseMethod("score")

###{{{ score.multigroup

score.multigroup <- function(x,data=x$data,p,...) {
  ## Check for random slopes
  Xfix <- FALSE
  xfix <- list()
  xx <- x
  for (i in 1:x$ngroup) {
    x0 <- x$lvm[[i]]
    data0 <- x$data[[i]]
    xfix0 <- colnames(data0)[(colnames(data0)%in%parlabels(x0))]
    xfix <- c(xfix, list(xfix0))
    if (length(xfix0)>0) { ## Yes, random slopes
      Xfix<-TRUE
    }
    xx$lvm[[i]] <- x0
  }   
  pp <- modelPar(xx,p)$p
  S <- 0; for (i in 1:x$ngroup)
    S <- S + colSums(score(x$lvm[[i]],p=pp[[i]],data=data[[i]]))
  return(S)  
}

###}}} score.multigroup

###{{{ score.lvmfit
score.lvmfit <- function(x, data=model.frame(x), p=pars(x), ...) {
  score(x$model0,data=data,p=p,...)
}

###}}} score.lvmfit

###{{{ score.lvm

score.lvm <- function(x, data, p, S, n, mu=NULL, weight=NULL, debug=FALSE, reindex=FALSE, mean=TRUE, constrain=FALSE,...) {
  if (is.null(index(x)$dA) | reindex)
    x <- updatelvm(x,zeroones=TRUE,deriv=TRUE)
##    index(x) <- reindex(x,zeroones=TRUE,deriv=TRUE) ## If not already done, calculate some relevant zero-one matrices and matrix derivatives

  if (!is.null(data)) {
    xfix <- colnames(data)[(colnames(data)%in%parlabels(x))]
    xconstrain <- intersect(unlist(lapply(constrain(x),function(z) attributes(z)$args)),manifest(x))
    
    if ((missing(S) |  nrow(data)<2 | length(xfix)>0 | length(xconstrain)>0 | !is.null(weight)) & constrain) {
      Debug("before deriv.", debug)

      Debug(xfix,debug)
      x0 <- x
      if (length(xfix)>0) ## Random slopes?
        {          
          Debug("random slopes...",debug)
          nrow <- length(vars(x))
          xpos <- lapply(xfix,function(y) which(regfix(x)$labels==y))
          colpos <- lapply(xpos, function(y) ceiling(y/nrow))
          rowpos <- lapply(xpos, function(y) (y-1)%%nrow+1)
          myfix <- list(var=xfix, col=colpos, row=rowpos)
          for (i in 1:length(myfix$var)) 
            for (j in 1:length(myfix$col[[i]])) 
              regfix(x0, from=vars(x0)[myfix$row[[i]][j]],to=vars(x0)[myfix$col[[i]][j]]) <-
                data[1,myfix$var[[i]]]
          ##rep(data[1,myfix$var[[i]]],length(myfix$row[[i]]))
          index(x0) <- reindex(x0,zeroones=TRUE,deriv=TRUE)
          yvars <- endogenous(x0)
        }
        pp <- modelPar(x0,p)        
        
        myfun <- function(ii) {
          if (length(xfix)>0)
            for (i in 1:length(myfix$var)) {
              ##            for (j in 1:length(myfix$col[[i]])) {
              index(x0)$A[cbind(myfix$row[[i]],myfix$col[[i]])] <- data[ii,myfix$var[[i]]]
              ##            }
            }
          ##          browser()
          return(score(x0,data=data[ii,,drop=FALSE], p=with(pp,c(meanpar,p)),weight=weight[ii,,drop=FALSE],constrain=FALSE) )
        }
      return(t(sapply(1:nrow(data),myfun)))
    }
    
### Standard SEM:
      pp <- modelPar(x,p)
    
      mp <- modelVar(x,p,data=data[1,])       
      C <- mp$C
      xi <- mp$xi

      iC <- Inverse(C,0,det=FALSE)
      D <- with(pp, deriv(x, meanpar=meanpar, p=p, mom=mp, mu=NULL)) ##, all=length(constrain(x))>0))
      ##      D <- with(pp, deriv(x, meanpar=meanpar, mom=mp, mu=NULL))
      Debug("after deriv.", debug)
      myvars <- (index(x)$manifest)
      if (NCOL(data)!=length(myvars)) {
        data <- subset(data,select=myvars)
      }

    score <- c()
    score0 <- -1/2*as.vector(iC)%*%D$dS
      for (i in 1:NROW(data)) {
        z <- as.numeric(data[i,])
        u <- z-xi
        if (!is.null(weight)) {
          W <- diag(as.numeric(weight[i,]))
          score <- rbind(score,
                         as.numeric(crossprod(u,iC%*%W)%*%D$dxi +
                                    -1/2*(as.vector((iC
                                                     - iC %*% tcrossprod(u)
                                                     %*% iC)%*%W)) %*% D$dS
                                    
                                    ##                                    -1/2*as.vector(iC%*%W)%*%D$dS +                                      
                                    ##                                    1/2*as.vector(iC %*% tcrossprod(u)
                                    ##                                                  %*% iC %*% W)%*%D$dS
                                    ))
        } else {
          score <- rbind(score,
                         as.numeric(score0 + crossprod(u,iC)%*%D$dxi +
                                    1/2*as.vector(iC%*%tcrossprod(u)%*%iC)%*%D$dS))
        }
##                        1/2*t(u)%*%(t(u)%*%iC %x% iC) %*% D$dS))
        ##         score <- rbind(score,
        ##                        score0 - 0.5 * (-2*crossprod(u,iC)%*%D$dxi -
        ##                                        t(u)%*%(t(u)%*%iC (tcrossprod(u)) iC) %*% D$dS))                             
      }; colnames(score) <- names(p)
      return(score)
  }
  
  
  mp <- modelVar(x,p)  

### Here the emperical mean and variance of the population are sufficient statistics:
  C <- mp$C
  xi <- mp$xi
  iC <- Inverse(C,0,det=FALSE)
  Debug("Sufficient stats.",debug)
  if (!is.null(mu)) {
    W <- tcrossprod(mu-xi)
    T <- S+W
  } else {
    T <- S
  }
  D <- deriv(x, meanpar=attributes(mp)$meanpar, mom=mp, p=p, mu=mu, mean=mean) ##, all=length(constrain(x))>0)
  
  vec.iC <- as.vector(iC)  
##  Grad <- n/2*t(D$dS)%*%as.vector((iC%*%T%*%iC)) -
##    n/2*t(crossprod(vec.iC,D$dS))
  Grad <- n/2*crossprod(D$dS, as.vector(iC%*%T%*%iC)-vec.iC)
  if (!is.null(mu)) # & mp$npar.mean>0)
    Grad <- Grad - n/2*crossprod(D$dT,vec.iC)
  return(as.numeric(Grad))
}

gaussian.score <- function(data,mu,S,dS,dmu) {
##  D <- deriv(x, meanpar=attributes(mp)$meanpar, mom=mp, mu=mu)  
}

###}}} score.lvm

