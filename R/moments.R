Moments <- function(x,p,data,conditional=TRUE,...) {
  
}

`moments` <-
  function(x,...) UseMethod("moments")

moments.lvmfit <- function(x, p=pars(x),...) moments(Model(x),p=p,...)

moments.lvm <- function(x, p, debug=FALSE, conditional=FALSE, data=NULL, ...) {
##  moments.lvm <- function(x, p, meanpar=NULL, conditional=FALSE, debug=FALSE,...) {
### p: model-parameters as obtained from e.g. ´startvalues`
###       i.e. vector of regression parameters and variance parameters
### meanpar: mean-parameters (optional)

  ii <- index(x)

  pp <- modelPar(x,p)
##  AP <- with(pp, matrices(x,p,meanpar=TRUE,...))
  AP <- with(pp, matrices(x,p,meanpar=meanpar,data=data,...))
#  A <- AP$A; P <- AP$P; v <- AP$v;
  P <- AP$P
  if (!is.null(AP$v)) {
    names(AP$v) <- ii$vars
  }
  ##  rownames(P) <- colnames(P) <- ii$vars
  npar <- ii$npar
  npar.reg <- ii$npar.reg

  J <- ii$J
  Jidx <- ii$obs.idx
  if (conditional) {
    ##    mynames <- index(x)$endo.idx
    J <- ii$Jy
    px <- ii$px 
    P <-  px%*% tcrossprod(P, px)
    Jidx <- ii$endo.idx
  } else {
    ##    mynames <- ii$vars
    ##    J <- ii$J ## Manifest variable selection matrix    
  }
  
  Im <- diag(nrow(AP$A))  
  if (ii$sparse) {
    IAi <- with(AP, as(solve(Im-t(A)),"sparseMatrix"))
##    IAi <- as(solve(Diagonal(nrow(A))-t(A)),"sparseMatrix")
    G <- as(J%*%IAi,"sparseMatrix")
  } else {
    IAi <- solve(Im-t(AP$A))
    G <- J%*%IAi
    ##G <- IAi[Jidx,,drop=FALSE]
  }
  
  xi <- NULL
  if (!is.null(AP$v)) {
    xi <- G%*%AP$v ## Model-specific mean vector
  }
  ##  rownames(xi) <- J%*%mynames

  Cfull <- as.matrix(IAi %*% tcrossprod(P,IAi))
  C <- as.matrix(J %*% tcrossprod(Cfull,J))
  ##  C <- Cfull[Jidx,Jidx,drop=FALSE]
##  rownames(C) <- colnames(C) <- mynames
  
  return(list(Cfull=Cfull, C=C, v=AP$v, xi=xi, A=AP$A, P=P, IAi=IAi, J=J, G=G, npar=npar, npar.reg=npar.reg, npar.mean=ii$npar.mean, parval=AP$parval, constrain.idx=AP$constrain.idx, constrainpar=AP$constrainpar))
}
