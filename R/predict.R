predict.lvmfit <- function(object,data=model.frame(object),p=pars(object),...) {
  predict(Model(object),p=p,data=data,mu=object$mu,S=object$S,...)
}

predict.lvm <- function(object,p,data,path=FALSE,...) {
  ## data = data.frame of exogenous variables
  if (!all(exogenous(object)%in%colnames(data))) stop("dataframe should contain exogenous variables")
  
  ##  pp <- modelPar(object,p)
  m <- moments(object,p)
  if (path) {
    Y <- endogenous(object,top=TRUE)
    X <- setdiff(manifest(object),Y)
  } else {
    X <- exogenous(object)
    Y <- setdiff(manifest(object), X)
  }
  X.idx <- match(X,manifest(object))
  eta.idx <- match(latent(object),vars(object))
  obs.idx <- match(manifest(object),vars(object))
  ##  exo.idx <- match(exogenous(object),vars(object))
  ##  muX <- mu[X.idx]
  ##  varX <- S[X.idx, X.idx]

  X.idx.all <- match(X, vars(object))
  Y.idx.all <- match(Y, vars(object))

  ## Calculation of conditional variance given X=x
  A <- t(m$A); 
  IAi <- solve(diag(nrow(A))-A)
  P.x <- m$P; P.x[X.idx.all, X.idx.all] <- 0
  C.x <- (IAi%*% P.x %*%t(IAi))
  Cy.x <- C.x[Y.idx.all,Y.idx.all,drop=FALSE]

  px <- diag(nrow(A)); px[X.idx.all,X.idx.all] <- 0 ## Sets exogenous-entries to zero
  ##  print(sum(index(object)$px-px))
  
  other.idx <- match(setdiff(vars(object),endogenous(object)),vars(object))
  Jy <- diag(nrow(A)); Jy[,other.idx] <- 0; Jy <- Jy[-other.idx,,drop=FALSE]## Only keep endogenous entries
  ##  print(sum(Jy-index(object)$Jy))
  
  Cy.x2 <- Jy%*% (IAi %*% (px%*%m$P%*%t(px)) %*% t(IAi)) %*% t(Jy)
  
  
  ## Calculation of conditional mean given X=x
  G <- m$J%*%m$IAi
  mu.0 <- m$v; mu.0[X.idx.all] <- 0
  xs <- data[,X,drop=FALSE]
  mu.x <- apply(xs, 1, FUN=function(i) {res <- rep(0,length(mu.0)); res[X.idx.all] <- i; res})
  xi.x <- (m$IAi%*%(mu.0 + mu.x))

  Ey.x <- xi.x[Y.idx.all,,drop=FALSE]
  Eeta.x <- xi.x[eta.idx,,drop=FALSE]

  Ey.x2 <- t( Jy %*% ( m$IAi%*% (as.numeric(px%*%m$v) + mu.x)) )
##  print(t(Ey.x) - Ey.x2)
##  return(0)



  if (length(eta.idx)>0) {
    Ceta.x <- C.x[eta.idx,eta.idx]
    Lambda <- matrix(A[Y.idx.all,eta.idx,drop=FALSE], ncol=length(eta.idx))
    Cetay.x <- Ceta.x%*%t(Lambda)

    ys <- data[,Y]

    KK <- Cetay.x %*% solve(Cy.x)
    Eeta.y <- Eeta.x + KK %*% (t(ys) - (Ey.x))##)Ey.x)
##    Eeta.y <- Eeta.x + KK %*% (Ey.x-)
    
    
    Ceta.y <- Ceta.x - KK%*% t(Cetay.x)
  } else {
    Eeta.y <- NA
    Ceta.y <- NA
  }
  
  Vhat <- matrix(0, nrow(data), length(vars(object))); colnames(Vhat) <- vars(object)
  Vhat[,obs.idx] <- as.matrix(data[,manifest(object)])
  if (length(eta.idx)>0)
    Vhat[,latent(object)] <- t(Eeta.y)
  ##dd[,exogenous(m1),drop=FALSE]
  I <- diag(nrow(A));
  epsilonhat <- t( (I-A)%*%t(Vhat) - m$v )

  mydata <- matrix(0,ncol=ncol(A),nrow=nrow(data)); colnames(mydata) <- vars(object)
  mydata[,manifest(object)] <- as.matrix(data[,manifest(object)])
  Yhat <- t(mydata%*%t(A)) + (m$v)
##  Yhat <- t((as.matrix(data[,manifest(object)]))%*%t(A)) + (m$v)
  res <- t(Yhat)
##  res <- Ey.x ## Conditional mean
  
  attr(res, "cond.var") <- Cy.x
  attr(res, "blup") <- t(Eeta.y)
  attr(res, "var.blup") <- Ceta.y
  attr(res, "Ey.x") <- Ey.x
  attr(res, "eta.x") <- Eeta.x
  attr(res, "epsilon.y") <- epsilonhat
##  return(list(var.blup=Ceta.y, blup=t(Eeta.y), cond.var=Cy.x, cond.mean=t(Ey.x)))
  return((res))
}
