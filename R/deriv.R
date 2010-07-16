
deriv.lvm <- function(x, p, mom, cond=FALSE, meanpar=TRUE, mu=NULL, S=NULL, second=FALSE, zeroones=FALSE, all=!missing(mom),...) {
  ## if (!missing(mom)) {
  ##   if (mom$npar==length(attributes(mom)$pars))
  ##     meanpar <- NULL
  ## }
  if (missing(mom) & !missing(p)) {
    mom <- modelVar(x,p)
  }
  if (missing(mom) & !missing(p)) {
    mom <- modelVar(x,p)
    all <- TRUE
    if (mom$npar==length(p))
      meanpar <- NULL  
  }

  ii <- index(x)
  npar.total <- npar <- ii$npar; npar.reg <- ii$npar.reg
  npar.mean <- ifelse(is.null(meanpar),0,ii$npar.mean)
  if (npar.mean>0) {
    meanpar <- 1:npar.mean
  } else {
    meanpar <- NULL
  }

  nn <- x$parpos
  if (is.null(nn))  
    {
      nn <- matrices(x,1:npar + npar.mean,meanpar);
      nn$A[ii$M0!=1] <- 0
      nn$P[ii$P0!=1] <- 0
      nn$v[ii$v0!=1] <- 0
    }
  
  if (npar.reg>0) {
    regr.idx <- 1:npar.reg + npar.mean
  } else {
    regr.idx <- NULL
  }
  if (npar>(npar.reg)) {
    var.idx <- 1:(npar-npar.reg) + (npar.mean + npar.reg) ##(npar.reg+1):npar
  } else {
    var.idx <- NULL
  }
  
  if (!is.null(meanpar)) {
    mean.idx <- 1:npar.mean
    npar.total <- npar+length(mean.idx)
  }
  
  if (zeroones | is.null(ii$dA)) {
    dimA <- length(nn$A)
    if (ii$sparse) { ## Not used yet...
      if (!require("Matrix")) stop("package Matrix not available")
      dP <- dA <- Matrix(0, nrow=dimA, ncol=npar.total)
    } else {
      dP <- dA <- matrix(0, nrow=dimA, ncol=npar.total)
    }
    if (npar.reg>0) {
      dA[,regr.idx] <- sapply(regr.idx, function(i) izero(which(t(nn$A)==i),nrow(dA)) )
    }   
    par.var <- (npar.reg+1):npar
    if (npar>npar.reg) {
      dP[,var.idx] <- sapply(var.idx, function(i) izero(which(nn$P==i),nrow(dA)) )
    }    
    res <- list(dA=dA, dP=dP)
    
    if (!is.null(meanpar) & npar.mean>0) {
      if (ii$sparse) {
        dv <- Matrix(0, nrow=length(x$mean), ncol=npar.total)
      } else {
        dv <- matrix(0, nrow=length(x$mean), ncol=npar.total)
      }
      dv[,mean.idx] <- sapply(1:npar.mean, function(i) izero(which(nn$v==i),length(x$mean)) ) 
      res <- c(res, list(dv=dv))
    }
  } else {
    res <- with(ii, list(dA=dA, dP=dP, dv=dv))
  }
  if (!all) return(res)
 
  ## Non-linear constraints:
  if (!missing(p)  && length(index(x)$constrain.par)>0) {
##    print("Constraints....")
    for (pp in index(x)$constrain.par) {
      myc <- constrain(x)[[pp]]
      parval <- mom$parval
      vals <- parval[attributes(myc)$args]
      fval <- myc(unlist(vals))
      if (!is.null(attributes(fval)$grad)) {
        Gr <- attributes(fval)$grad(unlist(vals))
      } else {
        if (!require("numDeriv")) stop("numDeriv or analytical derivatives needed!")
        Gr <- as.numeric(jacobian(myc, unlist(vals)))
      }
      mat.idx <- mom$constrain.idx[[pp]]
      for (jj in 1:length(vals)) {
        allpars <- c(nn$A[attributes(vals[[jj]])$reg.idx[1]],
                     nn$P[attributes(vals[[jj]])$cov.idx[1]],
                     nn$v[attributes(vals[[jj]])$m.idx[1]])
        if (!is.null(mat.idx$cov.idx))
          res$dP[mat.idx$cov.idx,allpars] <- Gr[jj]
        if (!is.null(mat.idx$reg.idx))
          res$dA[mat.idx$reg.tidx,allpars] <- Gr[jj]
        if (!is.null(res$dv) & !is.null(mat.idx$m.idx))
          res$dv[mat.idx$m.idx,allpars] <- Gr[jj]
      }     
    }
  }

  if (is.null(ii$Kkk)) {
    nobs <- nrow(ii$J)
    ii$Ik <- diag(nobs)
    ii$Im <- diag(ncol(ii$A))
    ##    ii$Kkk <- commutation(nobs,sparse=FALSE)
  }

  N <- NCOL(ii$A)
  K <- nobs
##  browser()
  if (N>15) {
    dG <- with(mom, matrix(0,prod(dim(G)),NCOL(res$dA)))
    for (i in 1:NCOL(dG)) { ## vec(ABC) = (C'xA)*vec(B)
      dG[,i] <- as.vector(with(mom, G%*%matrix(res$dA[,i],ncol=N)%*%IAi))
    }
    G1 <- G2 <- G3 <- matrix(0,K^2,NCOL(res$dA))
    tGP <- with(mom, t(G%*%P))
    for (i in 1:NCOL(G1)) {
      G1[,i] <- as.vector(matrix(dG[,i],ncol=NCOL(mom$G))%*%tGP)
      G3[,i] <- as.vector(with(mom, ((G)%*%matrix(res$dP[,i],ncol=NCOL(mom$P))%*%t(G))))
    }
    G2 <- G1[as.vector(matrix(1:(K^2),K,byrow=TRUE)),]
    dS <- G1+G2+G3
  } else {
    dG <- suppressMessages(with(mom, (t(IAi) %x% G) %*% (res$dA)))  
    MM <- suppressMessages(with(mom, (G%*%P %x% ii$Ik)))
    G1<- MM %*% (dG)
    ## Commutatation product K*X: 
    ##  G2 <- with(mom, ii$Kkk%*%(G1))
    G2 <- G1[as.vector(matrix(1:(K^2),K,byrow=TRUE)),]
    G3 <- with(mom, (G%x%G)%*%(res$dP))
    dS <- G1+G2+G3
  }
  res <- c(res, list(dG=dG, dS=dS))

##    if (!is.null(meanpar)) {
##  if (!is.null(mu)) {
  if (!is.null(mom$v)){
      dxi <-
        with(mom, (t(v)%x% ii$Ik)%*%dG)
      if (!is.null(res$dv))
        dxi <- dxi+ mom$G%*%res$dv
      res <- c(res, list(dxi=dxi))
      if (!is.null(mu)) {
        muv <- mu-mom$xi
        dT <- suppressMessages(-(ii$Ik%x%muv + muv%x%ii$Ik) %*% dxi)
        res <- c(res, list(dT=dT))
      }      
    }
  
    if (second) {
      cat("Bleeding edge...\n")
      k <- nrow(ii$A)
      K <- ii$Kkk ## commutation(k,k)
      I <- ii$Ik ## diag(k)
      I2 <- diag(k*k)
      d2S1 <-  t(
                (I %x% K %x% I) %*% (
                                     ( I2 %x% as.vector(mom$G) )%*% dG +
                                     ( as.vector(mom$P) %x% I2 )%*% (dP)
                                     ) %*% t(dG)
                 )
      d2S2 <- K%*%d2S1 ### HK?
      d2S3 <- t(
                (I %x% K %x% I) %*% (
                                     ( I2 %x% as.vector(mom$G) )%*% dG +
                                     ( as.vector(mom$G) %x% I2 )%*% dG
                                     ) %*% t(dP)                
                )
      vec.d2S <- d2S1+d2S3+d2S3      
      res <- c(res, list(d2vecS=vec.d2S))
    }
  
  return(res)
}

