`matrices` <-
  function(x,...) UseMethod("matrices")

###{{{ matrices.lvm

matrices.lvm <- function(x,pars,meanpar=NULL,data=NULL,...) {
  ii <- index(x)
  A <- ii$A ## Matrix with fixed parameters and ones where parameters are free
  J <- ii$J ## Manifest variable selection matrix
  M0 <- ii$M0 ## Index of free regression parameters
  M1 <- ii$M1 ## Index of free and _unique_ regression parameters
  P <- ii$P  ## Matrix with fixed variance parameters and ones where parameters are free
  P0 <- ii$P0 ## Index of free variance parameters
  P1 <- ii$P1 ## Index of free and _unique_ regression parameters

  P1.lower <- P1[lower.tri(P1)]
  npar.reg <- sum(M1)
  npar.var <- sum(diag(P1)) + sum(P1.lower) ## Number of covariance parameters
  npar.mean <- ii$npar.mean

  constrain.par <- names(constrain(x))
  parval <- list()
    
  ##parname <- parname.all <- c()
  parname.all <- unique(x$par[!is.na(x$par)])
  parname <- setdiff(parname.all,constrain.par)

  if (npar.reg>0) {
    A[which(M1==1)] <- pars[seq_len(npar.reg)]
    for (p in parname) {
      idx <- which((x$par==p))
      newval <- A[idx[1]]
      attributes(newval)$reg.idx <- idx
      attributes(newval)$reg.tidx <- which(t(x$par==p))
      parval[[p]] <- newval
      if (length(idx)>1) {
        A[idx[-1]] <- parval[[p]]
      }
    } ## duplicate parameters
  }

  ##    browser()

  if (npar.reg==0) {
    pars.var <- pars
  } else {
    pars.var <- pars[-c(1:npar.reg)]
  }
  which.diag <- diag(P1==1)
  diag(P)[which.diag] <- pars.var[seq_len(sum(which.diag))]
  covparname.all <- unique(x$covpar[!is.na(x$covpar)])
  covparname <- setdiff(covparname.all,constrain.par)

  pars.off.diag <- pars.var
  if (sum(which.diag)>0) {
    pars.off.diag <- pars.off.diag[-c(1:sum(which.diag))]
  }
  counter <- 0
  if (length(pars.off.diag)>0 & ncol(P)>1)
  for (i in 1:(ncol(P1)-1))
    for (j in (i+1):nrow(P1)) {
      if (index(x)$P1[j,i]!=0) {
        counter <- counter+1
        P[j,i] <- pars.off.diag[counter]
      }
    }

  if (length(covparname)>0)
  for (p in covparname) {
    idx <- which(x$covpar==p)
    if (!(p%in%parname)) {
      parval[[p]] <- P[idx[1]]      
    }
    attributes(parval[[p]])$cov.idx <- idx
    if (length(idx)>1)
      P[idx[-1]] <- parval[[p]]
    if (npar.reg>0)
      if (p%in%parname) {
        idx.reg <- which(x$par==p)
        P[idx] <- A[idx.reg[1]]

        atr <- attributes(parval[[p]])
        parval[[p]] <- A[idx.reg[1]] ###?????
        attributes(parval[[p]]) <- atr

      }
  } ## duplicate parameters
  P[upper.tri(P)] <- t(P)[upper.tri(P)]    
##  P <- symmetrize(P)

  
  v <- NULL
  mparname.all <- NULL
##  if (!(is.null(meanpar)) | npar.mean!=0) {
##  if (!is.null(meanpar) | npar.mean==0)
    {    
    named <- sapply(x$mean, function(y) is.character(y) & !is.na(y))    
    fixed <- sapply(x$mean, function(y) is.numeric(y) & !is.na(y))
    v <- rep(0,length(x$mean))
    names(v) <- colnames(P)
    if (!is.null(meanpar) | npar.mean==0)    
      v[index(x)$v1==1] <- meanpar
    if (any(fixed))
        v[fixed] <- unlist(x$mean[fixed])
    mparname.all <- unique(x$mean[named])
    mparname <- setdiff(mparname.all,constrain.par)
    for (p in mparname) {
      idx <- which(x$mean==p)

      if (!(p%in%c(parname,covparname))) {
        parval[[p]] <- v[idx[1]]        
      }
      attributes(parval[[p]])$m.idx <- idx
            
      if (length(idx)>1)
        v[idx[-1]] <- parval[[p]]
      if (p %in% covparname & !(p %in% parname)) {
        idx.2 <- which(x$covpar==p)
        v[idx] <- P[idx.2[1]]
        ##  atr <- attributes(parval[[p]])
        ##        parval[[p]] <- P[idx.2[1]] ###?????
        ##        attributes(parval[[p]]) <- atr
      }
      if (p %in% parname) {
        idx.2 <- which(x$par==p)
        v[idx] <- A[idx.2[1]]
        ##        atr <- attributes(parval[[p]])
        ##        parval[[p]] <- A[idx.2[1]] ###?????
        ##        attributes(parval[[p]]) <- atr
      }
    } 
  }

  ## Constrained...
  constrain.idx <- NULL
  cname <- constrainpar <- c()  
  if (length(constrain.par)>0 && is.numeric(c(pars,meanpar))) {   
    constrain.idx <- list()
    for (p in constrain.par) {
      cname <- c(cname,p)
      reg.tidx <- reg.idx <- cov.idx <- m.idx <- NULL
      myc <- constrain(x)[[p]]
      xargs <- manifest(x)[na.omit(match(attributes(myc)$args,manifest(x)))]
      if (length(xargs)>0) {
        if (!is.null(data)) {
          parval[xargs] <- (data)[xargs]
        } else parval[xargs] <- 0
      }
      val <- unlist(c(parval,constrainpar,x$mean)[attributes(myc)$args])      
      cpar <- myc(val); 
      constrainpar <- c(constrainpar,list(cpar)); names(constrainpar) <- cname
      if (p%in%parname.all) {
        reg.idx <- which(x$par==p)
        reg.tidx <- which(t(x$par==p))
        if (!is.null(val))
          A[reg.idx] <- cpar##myc(val)
      }
      if (p%in%covparname.all) {
        cov.idx <- which(x$covpar==p)
        if (!is.null(val))
          P[cov.idx] <- cpar##myc(val)        
      }
      ##      if (!is.null(meanpar))
      if (p%in%mparname.all) {
          m.idx <- which(x$mean==p)
          if (!is.null(val))
            v[m.idx] <- cpar##myc(val)        
        }
      constrain.idx[[p]] <- list(reg.idx=reg.idx,reg.tidx=reg.tidx,cov.idx=cov.idx,m.idx=m.idx)
    }
  }
  
  if (x$index$sparse & !is.character(class(pars)[1])) {
    A <- as(A,"sparseMatrix")
    P <- as(P,"sparseMatrix")
    v <- as(v,"sparseMatrix")
  }
  return(list(A=A, P=P, v=v, parval=parval, constrain.idx=constrain.idx, constrainpar=constrainpar))
}

###}}} matrices.lvm

###{{{ matrices.multigroup

matrices.multigroup <- function(x, p) {
  pp <- modelPar(x,p)
  res <- list()
  for (i in seq_len(x$ngroup))
    res <- c(res, list(matrices(x$lvm[[i]],p=pp$par[[i]],meanpar=pp$mean[[i]])))
  return(res)
}

###}}}
