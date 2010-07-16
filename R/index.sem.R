parpos <- function(x,mean=TRUE,...) {
  if (mean)
    nn <- with(index(x),matrices(x,1:npar+npar.mean,meanpar=1:npar.mean)) ## Position of parameters
  else nn <- with(index(x),matrices(x,1:npar,NULL))
  nn$A[index(x)$M0!=1] <- 0
  nn$P[index(x)$P0!=1] <- 0
  nn$v[index(x)$v0!=1] <- 0
  nn
}

updatelvm <- function(x,mean=TRUE,...) {
  index(x) <- reindex(x,mean=mean,...)
  x$parpos <- parpos(x,mean=mean,...)
  return(x)
}

"index" <- function(x,...) UseMethod("index")
"index<-" <- function(x,...,value) UseMethod("index<-")

"index.lvm" <- function(x,...) { x$index }
"index.lvmfit" <- function(x,...) { index(Model(x)) }

"index<-.lvm" <- function(x,...,value)  { x$index <- value; return(x) }
"index<-.lvmfit" <- function(x,...,value) { Model(x)$index <- value; return(x) }


###   A  ## Matrix with fixed parameters and ones where parameters are free
###   J  ## Manifest variable selection matrix
###   M0 ## Index of free regression parameters
###   M1 ## Index of free and _unique_ regression parameters
###   P  ## Matrix with fixed variance parameters and ones where parameters are free
###   P0 ## Index of free variance parameters
###   P1 ## Index of free and _unique_ regression parameters
###   npar.var  ## Number of covariance parameters
`reindex` <-
function(x, debug=FALSE,sparse=FALSE,standard=TRUE,zeroones=FALSE,deriv=FALSE,mean=TRUE) { ## Extract indices of parameters from model
  M <- as(Graph(x), Class="matrix")
  Debug("M=",debug)

  eta <- latent(x) ## Latent variables/Factors
  m <- length(eta)
  obs <- manifest(x)  ## Manifest/Observed variables
  endo <- endogenous(x)
  exo <- exogenous(x)
  allvars <- vars(x)
  eta.idx <- match(eta,allvars)
  obs.idx <- match(obs,allvars)
  exo.idx <- match(exo,allvars)
  exo.obsidx <- match(exo,obs)
  endo.obsidx <- match(endo,obs)

  fix.idx <- !is.na(x$fix) ## Index of fixed parameters
  covfix.idx <- !is.na(x$covfix) ## Index of fixed covariance parameters

  constrain.par <- NULL
  if (length(constrain(x))>0) constrain.par <- names(constrain(x))

  Debug("M0",debug)
  M0 <- M;  M0[fix.idx] <- 0 ## Matrix of indicators of free regression-parameters (removing fixed parameters)
  M1 <- M0; ## Matrix of indiciator of free _unique_ regression parameters (removing fixed _and_ duplicate parameters)
  parname <- unique(x$par[!is.na(x$par)])
##  parname.all <- unique(x$par[!is.na(x$par)])
##  parname <- setdiff(parname.all,constrain.par)
  for (p in parname) {
    ii <- which(x$par==p)
    if (length(ii)>1)
      M1[ii[-1]] <- 0
    if (p %in% constrain.par)
      M0[ii] <- M1[ii] <- 0
  }
  ## if (length(constrain(x))>0) {
  ##   jj <- which(x$par%in%constrain.par)
  ##   M1[jj] <- 0
  ## }  
  npar.reg <- sum(M1) ## Number of free regression parameters
  Debug("npar done",debug)

  P <- x$cov;
  P0 <- P;  P0[covfix.idx] <- 0 ## Matrix of indicators of free covariance-parameters (removing fixed parameters)  
  P1 <- P0 ## Matrix of indiciator of free _unique_ variance parameters (removing fixed _and_ duplicate parameters)
  covparname <- unique(x$covpar[!is.na(x$covpar)])
##  covparname.all <- unique(x$covpar[!is.na(x$covpar)])
##  covparname <- setdiff(covparname.all,constrain.par)
  for (p in covparname) {
    ii <- which(x$covpar==p)
    if (length(ii)>1)
      P1[ii[-1]] <- 0
    if (p%in%c(parname,constrain.par))
      P0[ii] <- P1[ii] <- 0    
  } 
  ##  P1[upper.tri(P1)] <- 0
  ##  P1 <- symmetrize(P1) ### OBS CHECK ME
  
  npar.var <- sum(c(diag(P1),P1[lower.tri(P1)]))
  parnames <- paste("p", 1:(npar.reg+npar.var), sep="")
  
  Debug(npar.reg,debug)
  
##  A <- M0;
  A <- M
  A[fix.idx] <- x$fix[fix.idx] ## ... with fixed parameters in plac
##  P <- P0;
  P[covfix.idx] <- x$covfix[covfix.idx] ## ... with fixed parameters in plac

  px <- Jy <- J <- I <- diag(length(vars(x)))
  if (m>0) {
    J[eta.idx,eta.idx] <- 0; J <- J[-eta.idx,,drop=FALSE]
  } ## Selection matrix (selecting observed variables)
  {
    ## Selection matrix (selection endogenous variables)
    Jy[c(eta.idx,exo.idx),c(eta.idx,exo.idx)] <- 0; Jy <- Jy[-c(eta.idx,exo.idx),,drop=FALSE]
    ## Cancelation matrix (cancels rows with exogenous variables)
    px[exo.idx,exo.idx] <- 0
  } 

  ## Creating indicitor of free mean-parameters
  fixed <- sapply(x$mean, function(y) is.numeric(y) & !is.na(y))
  named <- sapply(x$mean, function(y) is.character(y) & !is.na(y))
  mparname <- unique(x$mean[named])
##  mparname.all <- unique(x$mean[named])
##  mparname <- setdiff(mparname.all,constrain.par)
  v0 <- rep(1,length(x$mean)) ## Vector of indicators of free mean-parameters
  v0[fixed] <- 0; v1 <- v0
  for (p in mparname) {
    idx <- which(x$mean==p)
    if (length(idx)>1) {
##      print(idx[-1])
      v1[idx[-1]] <- 0
    }
    if (p%in%c(parname,covparname,constrain.par))
      v0[idx] <- v1[idx] <- 0
  } ## duplicate parameters
  
  
  ## Return:
  ## Adjacency-matrix (M)
  ## Matrix of regression-parameters (0,1) _with_ fixed parameters (A)
  ## Matrix of variance-parameters (indicators 0,1) (P)
  ## Manifest selection matrix (J),
  ## Position of variables matrix (Apos),
  ## Position of covariance variables matrix (Ppos),
  ## Position/Indicator matrix of free regression parameters (M0)
  res <- list(vars=allvars, manifest=obs, exogenous=exo, latent=eta,
              endo=endo,
              exo.idx=exo.idx, eta.idx=eta.idx,
              exo.obsidx=exo.obsidx, endo.obsidx=endo.obsidx,
              obs.idx=obs.idx, 
              endo.idx=setdiff(obs.idx,exo.idx))
  
  if (standard) {
    res <- c(res, list(M=M, A=A, P=P, P0=P0, P1=P1, M0=M0, M1=M1, v0=v0, v1=v1, npar=(npar.reg+npar.var), npar.reg=npar.reg, npar.mean=sum(v1), constrain.par=constrain.par))
    npar.total <- res$npar+res$npar.mean
    ##  if (sparse)
    ##    res <- lapply(res, function(x) if(is.matrix(x)) as(x,"sparseMatrix") else x)
    res <- c(res, list(J=J, Jy=Jy, px=px, sparse=sparse))
  } else {
    res <- index(x)
  }
  
  if (zeroones) {
    if (sparse) {
      if (!require("Matrix")) stop("package Matrix not available")
      Ik <- Matrix:::Diagonal(length(obs))
      Im <- Matrix:::Diagonal(ncol(A))
      ##      Kkk <- commutation(length(obs),sparse=TRUE)
      Kkk <- NULL
      J <- as(J, "sparseMatrix")
      Jy <- as(Jy, "sparseMatrix")
      px <- as(px, "sparseMatrix")
      
    } else {
      Ik <- diag(length(obs))
      Im <- diag(ncol(A))
      ##      Kkk <- commutation(length(obs),sparse=FALSE)
    }
    Kkk <- NULL

    
    res[c("Ik","Im","Kkk")] <- NULL
    res <- c(res, list(Ik=Ik, Im=Im, Kkk=Kkk))
  }
  if (deriv) {
##    print("DERIV")
    if (res$npar.mean>0 & mean)
      D <- deriv(x,meanpar=rep(1,res$npar.mean),zeroones=TRUE)
    else
      D <- deriv(x,meanpar=NULL,zeroones=TRUE)
    res[c("dA","dP","dv")] <- NULL
    res <- c(res, list(dA=D$dA, dP=D$dP, dv=D$dv))
  }
  Debug("index done",debug)
  
  return(res)
}

