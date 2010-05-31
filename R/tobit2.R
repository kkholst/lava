##  Estimator for right-censored data

tobit2_objective.lvm <- function(x,p,data,weight,indiv=FALSE,debug=FALSE,algorithm=GenzBretz(abseps=1e-5),...) {
  require(mvtnorm)
  zz <- manifest(x)
  Debug("start",debug)
  d <- as.matrix(data[,zz,drop=FALSE]);
  colnames(d) <- zz
  Debug("dd",debug)
  yy <- endogenous(x)
  yy.idx <- match(yy,zz)
  Status <- matrix(1,ncol=length(zz),nrow=nrow(d))
  Status[,yy.idx]<- as.matrix(weight)[,yy,drop=FALSE]

  Debug("dd",debug)

  patterns <- unique(Status,MARGIN=1)
  cens.type <- apply(Status,1,
                     function(x) which(apply(patterns,1,function(y) identical(x,y))))  
  mp <- modelVar(x,p)
  Sigma <- mp$C ## Model specific covariance matrix
  xi <- mp$xi ## Model specific mean-vector
  Debug("loop:",debug)
  
##  val <- 0
  val <- c()
  for (i in 1:nrow(patterns)) {
    ## Usual marginal likelihood for status==1
    pat <- patterns[i,]
    idx <- which(cens.type==i)
    noncens.idx <- which(pat==1)
    cens.idx <- which(pat==0)
    noncens.y <- zz[noncens.idx]
    cens.y <- setdiff(zz,noncens.y)
    y <- d[idx,,drop=FALSE]
    val1 <- 0; val0 <- 0
    if (length(noncens.y)>0) {
      ## p(y), using: int[p(y,y*)]dy* =  p(y) int[p(y*|y)]dy*
      val1 <- dmvnorm(d[idx,noncens.y,drop=FALSE], ##-mu[idx,noncens.idx,drop=FALSE],
                      mean=xi[noncens.idx],
                      sigma=Sigma[noncens.idx,noncens.idx,drop=FALSE],log=TRUE)
    }
    if (length(cens.y)>0) {
      if (length(noncens.y)==0) {
        val0 <- sapply(idx,
                       function(ii) log(pmvnorm(lower=as.numeric(d[ii,cens.idx]),
                                                mean=as.numeric(xi),
                                                sigma=Sigma,algorithm=algorithm)))
      } else {
         M <- mom.cens(x,p,data=y,cens.idx=cens.idx,conditional=TRUE,deriv=FALSE)
         val0 <- c()
         for (j in 1:length(idx)) {
           val0 <- c(val0,log(pmvnorm(lower=as.numeric(y[j,cens.idx]),mean=M$mu.censIobs[j,],sigma=M$S.censIobs,algorithm=algorithm)) )
         }        
        }
    }
    val <- c(val,-val1-val0)
  }
  if (!indiv)
    return(sum(val))
  val
}
tobit2_logLik.lvm <- function(object,p,data,weight, ...) {
  res <- -tobit2_objective.lvm(x=object,p=p,data=data,weight=weight,...) - gaussian_logLik.lvm(object,p=p,data=data,type="exo",weight=NULL,...)
  return(res)
}

cens.score2 <- function(x,p,data,cens.idx,cens.which.left,...) {
  obs.idx <- setdiff(1:NCOL(data),cens.idx)
  n <- NROW(data)
  M <- mom.cens(x,p,data=data,cens.idx=cens.idx,conditional=TRUE,deriv=TRUE)

#  L <- diag(length(cens.idx))
#  L[cens.which.left,cens.which.left] <- (-1)  
##  print(M)
  ## Censored part:
  if (length(cens.idx)>0) {
    A <- (-1)
##     mu <- M$mu.cens
##     S <- M$S.cens
##     dmu <- M$dmu.cens
##     dS <- M$dS.cens
    z <- A*matrix(data[,cens.idx],nrow=n)
    ##    z <- data[,cens.idx,drop=FALSE]
    DCDFs <- c()
    for (i in 1:n) {
      DCDF <- Dthetapmvnorm(z[i,],
                            mu=A*M$mu.censIobs[i,],
                            S=M$S.censIobs,
                            dS=M$dS.censIobs,
                            dmu=A*matrix(M$dmu.censIobs[,,i],nrow=length(cens.idx)),
                            ...)
      alpha <- attributes(DCDF)$cdf
      DCDFs <- rbind(DCDFs, 1/alpha*DCDF)
    }
    ##    DCDF <- Dthetapmvnorm(z,mu=M$mu.cens,S=S,dmu=M$dmu.cens,dS=M$dS.cens)
    ##    S0 <- 1/attributes(DCDF)$cdf*DCDF
    S0 <- DCDFs    
  } else S0 <- 0
  ## Observed part:
  if (length(obs.idx)>0) {
    y1 <- as.matrix(data[,obs.idx,drop=FALSE])
    if (n==1) y. <- y1 else y. <- colMeans(y1)
    mu <- M$mu.obs
    S <- M$S.obs
    iS <- solve(S)
    dS <- M$dS.obs
    dmu <- M$dmu.obs
    S1 <- c()
    part0 <- -1/2*as.vector(iS)%*%dS
    for (i in 1:n) {
      z <- as.numeric(y1[i,])
      u <- z-mu
      S1 <- rbind(S1,
                  as.numeric(part0 + crossprod(u,iS)%*%dmu +
                             1/2*as.vector(iS%*%tcrossprod(u)%*%iS)%*%dS))
    }
  } else S1 <- 0
##   cat(rep("*",50),sep="")
##   print(S1)
##   cat(rep("-",50),sep="")
##   print(S0)
##   cat(rep("#",50),sep="")
##   print(dim(data))
  return((S1+S0))
}


tobit2_hessian.lvm <- function(p,model,data,weight,...) {
  S <- -tobit2_gradient.lvm(model,p=p,data=data,weight=weight,indiv=TRUE,...)
  J <- t(S)%*%S
  attributes(J)$grad <- colSums(S)
  return(J)  
}
tobit2_gradient.lvm <- function(x,p,data,weight,indiv=FALSE,...) {
  require(mvtnorm)
  zz <- manifest(x)
  d <- as.matrix(data[,zz,drop=FALSE]); colnames(d) <- zz
  yy <- endogenous(x)
  yy.idx <- match(yy,zz)
  Status <- matrix(1,ncol=length(zz),nrow=nrow(d))
  Status[,yy.idx] <- as.matrix(weight)[,yy,drop=FALSE]
  patterns <- unique(Status,MARGIN=1)
  cens.type <- apply(Status,1,
                     function(x) which(apply(patterns,1,function(y) identical(x,y))))  
  val <- 0
  score <- c()
  for (i in 1:nrow(patterns)) {
    ## Usual marginal likelihood for status==1
    pat <- patterns[i,]
    idx <- which(cens.type==i)
    noncens.idx <- which(pat==1)
    cens.idx <- which(pat==0)
    y <- d[idx,,drop=FALSE]
    dummy <- cens.score2(x,p,data=y,cens.idx=cens.idx)
    score <- rbind(score,dummy)
  }
  if (indiv)
    return(-score)
  return(-colSums(score))
}

