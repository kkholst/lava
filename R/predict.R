##' @S3method predict lvmfit
predict.lvmfit <- function(object,x=NULL,data=model.frame(object),p=pars(object),...) { 
  predict(Model(object),x=x,p=p,data=data,...)
}

##' @S3method predict lvm
##' @aliases predict.lvmfit
##' @param object Model object
##' @param x optional list of (endogenous) variables to condition on
##' @param residual If true the residuals are predicted
##' @param p Parameter vector 
##' @param data Data to use in prediction
##' @param path Path prediction (...)
##' @param quick If TRUE the conditional mean and variance given covariates are returned (and all other calculations skipped)
##' @param ... Additional arguments to lower level function
##' @examples
##' m <- lvm()
##' @export 
predict.lvm <- function(object,x=NULL,residual=FALSE,p,data,path=FALSE,quick=is.null(x)&!(residual|path),...) {
  ## data = data.frame of exogenous variables

  if (!quick && !all(exogenous(object)%in%colnames(data))) stop("dataframe should contain exogenous variables")
  m <- moments(object,p,data=data)
  if (quick) { ## Only conditional moments given covariates
    ii <- index(object)
    P.x <- m$P; P.x[ii$exo.idx, ii$exo.idx] <- 0
    Cy.x <- (m$IAi%*% tcrossprod(P.x,m$IAi))[ii$endo.idx,ii$endo.idx,drop=FALSE]
    ## Cy.x2 <- m$C[ii$endo.obsidx,ii$endo.obsidx]-m$C[ii$endo.obsidx,ii$exo.obsidx]%*%solve(m$C[ii$exo.obsidx,ii$exo.obsidx])%*%m$C[ii$exo.obsidx,ii$endo.obsidx]
    X <- ii$exogenous
    mu.0 <- m$v; mu.0[ii$exo.idx] <- 0
    if (length(X)>0) {
      mu.x <- matrix(0,ncol=nrow(data),nrow=length(mu.0))
      mu.x[ii$exo.idx,] <- t(data[,X,drop=FALSE])
      xi.x <- (m$IAi[ii$endo.obsidx,]%*%(mu.0 + mu.x))
    } else {
      xi.x <- matrix(as.vector(m$IAi[ii$endo.obsidx,]%*%mu.0),ncol=nrow(data),nrow=length(mu.0))
      rownames(xi.x) <- names(mu.0)
    }
    return(structure(t(xi.x),cond.var=Cy.x,
                     p=m$p,
                     e=m$e))
  }
  
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
  X.idx.all <- match(X, vars(object))
  Y.idx.all <- match(Y, vars(object))

  ## Calculation of conditional variance given X=x
  A <- t(m$A);
  P <- m$P
  P.x <- m$P; P.x[X.idx.all, X.idx.all] <- 0
  C.x <- (m$IAi%*% P.x %*%t(m$IAi))
  Cy.x <- C.x[Y.idx.all,Y.idx.all,drop=FALSE]
  ##  px <- diag(nrow(A)); px[X.idx.all,X.idx.all] <- 0 ## Sets exogenous-entries to zero
  ##  print(sum(index(object)$px-px))
  ##  other.idx <- match(setdiff(vars(object),endogenous(object)),vars(object)) ##endogenous entries
  ## Calculation of conditional mean given X=x
  G <- m$J%*%m$IAi
  mu.0 <- m$v; mu.0[X.idx.all] <- 0
  if (length(X)>0) {
    xs <- data[,X,drop=FALSE]
    mu.x <- apply(xs, 1, FUN=function(i) {res <- rep(0,length(mu.0)); res[X.idx.all] <- i; res})
    xi.x <- (m$IAi%*%(mu.0 + mu.x))
  } else {
    xi.x <- matrix(as.vector(m$IAi%*%mu.0),ncol=nrow(data),nrow=length(mu.0))
    rownames(xi.x) <- names(mu.0)
  }
  Ey.x <- xi.x[Y.idx.all,,drop=FALSE]
  Eeta.x <- xi.x[eta.idx,,drop=FALSE]
  Cy.epsilon <- P.x%*%t(m$IAi)
  Czeta.y <- Cy.epsilon[eta.idx,endogenous(object)]
  ##  Eeta.x + t([,c(4,8)])[,endogenous(e)]%*%solve(Cy.x)%*%t(rr)
  
  ## m <- moments(object,p)
  ## S <- m$Cfull
  ## v <- m$v
  IA <- diag(nrow=nrow(m$A))-t(m$A)
  ## mu <- as.vector(m$IAi%*%m$v); names(mu) <- names(m$v)
  ##  S <- C.x
  ##  mu <- t(xi.x)

  ys <- data[,Y]
  ry <- t(ys)-Ey.x
  y <- NULL
  if (!is.null(x)) {
    if (class(x)[1]=="formula")  {
      xy <- getoutcome(x)
      if (length(xy)>0) {
        y <- decomp.specials(xy)
      }
      x <- attributes(xy)$x
    }
    if (length(x)==0) {
      if (!is.null(y))
        xi.x <- xi.x[y,,drop=FALSE]
      return(t(xi.x))
    }
    x <- intersect(x,endogenous(object))
    if (is.null(y))
      y <- setdiff(vars(object),c(x,exogenous(object)))

    E.x <- xi.x[y,] + C.x[y,x]%*%solve(C.x[x,x])%*%ry[x,]
    if (residual) {
      Vhat <- matrix(0, nrow(data), length(vars(object))); colnames(Vhat) <- vars(object)
      Vhat[,obs.idx] <- as.matrix(data[,manifest(object)])
      Vhat[,y] <- t(E.x)
      return(t((IA%*%t(Vhat)-m$v)))
    }
    res <- t(E.x); colnames(res) <- y
    return(res)
##      return(t(Czeta.y[,x]%*%solve(Cy.x[x,x])%*%ry[x,]))
##    return(t(Eeta.x + C.x[eta.idx,x]%*%solve(Cy.x[x,x])%*%ry[x,]))
  }

  if (length(eta.idx)>0) {
    Ceta.x <- C.x[eta.idx,eta.idx]
    Lambda <- matrix(A[Y.idx.all,eta.idx,drop=FALSE], ncol=length(eta.idx))
    Cetay.x <- Ceta.x%*%t(Lambda)
    KK <- Cetay.x %*% solve(Cy.x)
    Eeta.y <- Eeta.x + KK %*% ry ##(t(ys) - (Ey.x))##)Ey.x)
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
  I <- diag(nrow=nrow(A));
  epsilonhat <- (t( IA%*%t(Vhat) - m$v ))[,c(endogenous(object),latent(object)),drop=FALSE]
  if (residual) {
    return(epsilonhat)
  }
  
  mydata <- matrix(0,ncol=ncol(A),nrow=nrow(data)); colnames(mydata) <- vars(object)
  mydata[,manifest(object)] <- as.matrix(data[,manifest(object)])
  for (i in latent(object))
    mydata[,i] <- m$v[i]
  
  Yhat <- t(mydata%*%t(A)) + (m$v)
##  Yhat <- t((as.matrix(data[,manifest(object)]))%*%t(A)) + (m$v)
##  res <- t(Yhat)
  res <- t(Ey.x) ## Conditional mean
  
  ##  attr(res, "cond.var") <- t(Yhat)
  attr(res, "cond.var") <- Cy.x
  attr(res, "blup") <- t(Eeta.y)
  attr(res, "var.blup") <- Ceta.y
  attr(res, "Ey.x") <- Ey.x
  attr(res, "eta.x") <- Eeta.x
  attr(res, "epsilon.y") <- epsilonhat
  attr(res, "p") <- m$p
  attr(res, "e") <- m$e
##  return(list(var.blup=Ceta.y, blup=t(Eeta.y), cond.var=Cy.x, cond.mean=t(Ey.x)))
  class(res) <- c("lvm.predict","matrix")
  return(res)
}

##' @S3method print lvm.predict
print.lvm.predict <- function(x,...) print(x[,])
