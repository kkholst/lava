##' @S3method deriv lvm
deriv.lvm <- function(expr, p, mom, conditional=FALSE, meanpar=TRUE, mu=NULL, S=NULL, second=FALSE, zeroones=FALSE, all=!missing(mom),...) {

  if (missing(mom) & !missing(p)) {
    mom <- modelVar(expr,p,conditional=conditional,...)
    all <- TRUE
    if (mom$npar==length(p))
      meanpar <- NULL  
  }

  ii <- index(expr)
  npar.total <- npar <- ii$npar; npar.reg <- ii$npar.reg
  npar.mean <- ifelse(is.null(meanpar),0,ii$npar.mean)
  meanpar <- seq_len(npar.mean)
  
  nn <- expr$parpos
  if (is.null(nn))  
    {
      nn <- matrices.lvm(expr,1:npar + npar.mean,meanpar);
      nn$A[ii$M0!=1] <- 0
      nn$P[ii$P0!=1] <- 0
      nn$v[ii$v0!=1] <- 0
    }
 
  regr.idx <- seq_len(npar.reg) + npar.mean
  var.idx <- seq_len(npar-npar.reg) + (npar.mean + npar.reg)   
  mean.idx <- seq_len(npar.mean)
  npar.total <- npar+length(mean.idx)

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
    par.var <- seq(npar.reg+1,npar)
    if (npar>npar.reg) {
      dP[,var.idx] <- sapply(var.idx, function(i) izero(which(nn$P==i),nrow(dA)) )
    }    
    res <- list(dA=dA, dP=dP)
    
    {
      if (ii$sparse) {
        dv <- Matrix(0, nrow=length(expr$mean), ncol=npar.total)
      } else {
        dv <- matrix(0, nrow=length(expr$mean), ncol=npar.total)
      }
      if (!is.null(meanpar) & npar.mean>0)
        dv[,mean.idx] <- sapply(1:npar.mean, function(i) izero(which(nn$v==i),length(expr$mean)) ) 
      res <- c(res, list(dv=dv))
    }
  } else {
    res <- with(ii, list(dA=dA, dP=dP, dv=dv))
    for (pp in nn$parval) {
      res$dP[attributes(pp)$cov.idx,pp] <- 1
      res$dv[attributes(pp)$m.idx,pp] <- 1      
    }
    nn$parval
  }

  if (!all) return(res)
  ## Non-linear constraints:
  cname <- constrainpar <- c()
  if (!missing(p)  && length(index(expr)$constrain.par)>0) {
    for (pp in index(expr)$constrain.par) {
      myc <- constrain(expr)[[pp]]
      if (!is.null(myc)) {
        parval <- mom$parval
        vals <- c(parval,constrainpar,mom$v)[attributes(myc)$args]
           
        fval <- myc(unlist(vals))
        if (!is.null(attributes(fval)$grad)) {
          Gr <- attributes(fval)$grad(unlist(vals))
        } else {
          if (!require("numDeriv")) stop("numDeriv or analytical derivatives needed!")
          Gr <- as.numeric(jacobian(myc, unlist(vals)))
      }
        mat.idx <- mom$constrain.idx[[pp]]
        cname <- c(cname,pp)
        attributes(fval)$grad <- Gr
        attributes(fval)$vals <- vals
        constrainpar <- c(constrainpar,list(fval)); names(constrainpar) <- cname
        
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
  }
  if (is.null(ii$Kkk)) {
    nobs <- nrow(mom$J)
    ii$Ik <- diag(nobs)
    ii$Im <- diag(ncol(ii$A))
    ##    ii$Kkk <- commutation(nobs,sparse=FALSE)
  }

  N <- NCOL(ii$A)
  K <- nobs
  ## if (N>10) {
  dG <- with(mom, kronprod(t(IAi),G,res$dA))
  GP <- with(mom,G%*%P)
  G1 <- with(mom, kronprod(GP,ii$Ik,dG))
  G2 <- G1[as.vector(matrix(1:(K^2),K,byrow=TRUE)),]
  G3 <- with(mom, kronprod(G,G,res$dP))
  dS <- G1+G2+G3
  ## } else {
  ##   dG <- suppressMessages(with(mom, (t(IAi) %x% G) %*% (res$dA)))  
  ##   MM <- suppressMessages(with(mom, (G%*%P %x% ii$Ik)))
  ##   G1<- MM %*% (dG)
  ##   ## Commutatation product K*X: 
  ##   ##  G2 <- with(mom, ii$Kkk%*%(G1))
  ##   G2 <- G1[as.vector(matrix(1:(K^2),K,byrow=TRUE)),]
  ##   G3 <- with(mom, (G%x%G)%*%(res$dP))
  ##   dS <- G1+G2+G3
  ## }
  ## }
  res <- c(res, list(dG=dG, dS=dS))
  
  if (!is.null(mom$v)) {
      ## dxi <-        
      ##   with(mom, (t(v)%x% ii$Ik)%*%dG)
      dxi <-        
        with(mom, kronprod(t(v),ii$Ik,dG))
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
      k <- nrow(ii$A)
      K <- ii$Kkk ## commutation(k,k)
      I <- ii$Ik ## diag(k)
      I2 <- diag(k*k)
      ##      KI <- I[as.vector(matrix(1:(K^2),K,byrow=TRUE)),]
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

