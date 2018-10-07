##' @export
deriv.function <- function(expr, parameter_, ..., parameter_.increment =.Machine$double.xmin) {
    p <- length(parameter_)
    f0 <- expr(parameter_)
    z0 <- numeric(p)
    res <- matrix(NA,nrow=length(f0),ncol=p)
    for (i in seq(p)) {
        z <- z0; z[i] <- parameter_.increment*1i
        res[,i] <- Im(expr(parameter_+z,...))/parameter_.increment
    }
    res       
}

##' @export
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
    npar.ex <- ii$npar.ex
    meanpar <- seq_len(npar.mean)
    nn <- expr$parpos

    if (is.null(nn))
    {
        nn <- matrices2(expr, seq_len(npar+npar.mean+npar.ex));
        nn$A[ii$M0!=1] <- 0
        nn$P[ii$P0!=1] <- 0
        nn$v[ii$v0!=1] <- 0
        nn$e[ii$e0!=1] <- 0
    }

    regr.idx <- seq_len(npar.reg) + npar.mean
    var.idx <- seq_len(npar-npar.reg) + (npar.mean + npar.reg)
    mean.idx <- seq_len(npar.mean)
    npar.total <- npar+length(mean.idx)
    epar.idx <- seq_len(npar.ex)+npar.total
    npar.total <- npar.total+length(epar.idx)

    if (zeroones | is.null(ii$dA)) {
        dimA <- length(ii$A)
        if (ii$sparse) { ## Not used yet...
            if (!requireNamespace("Matrix",quietly=TRUE)) stop("package Matrix not available")
            dP <- dA <- Matrix::Matrix(0, nrow=dimA, ncol=npar.total)
        } else {
            dP <- dA <- matrix(0, nrow=dimA, ncol=npar.total)
        }
        if (npar.reg>0) {
            dA[,regr.idx] <- sapply(regr.idx, function(i) izero(which(t(nn$A)==i),nrow(dA)) )
        }
        if (npar>npar.reg) {
            dP[,var.idx] <- sapply(var.idx, function(i) izero(which(nn$P==i),nrow(dA)) )
        }
        res <- list(dA=dA, dP=dP)

        {
            if (ii$sparse) {
                dv <- Matrix::Matrix(0, nrow=length(expr$mean), ncol=npar.total)
            } else {
                dv <- matrix(0, nrow=length(expr$mean), ncol=npar.total)
            }
            if (!is.null(meanpar) & npar.mean>0)
                dv[,mean.idx] <- sapply(mean.idx, function(i) izero(which(nn$v==i),length(expr$mean)) )
            res <- c(res, list(dv=dv))
        }
    } else {
        res <- with(ii, list(dA=dA, dP=dP, dv=dv))
        for (pp in nn$parval) {
            res$dP[attributes(pp)$cov.idx,pp] <- 1
            res$dv[attributes(pp)$m.idx,pp] <- 1
        }
    }

    if (!all) return(res)
    ## Non-linear constraints:
    cname <- constrainpar <- c()
    if (!missing(p) && length(index(expr)$constrain.par)>0) {
        for (pp in index(expr)$constrain.par) {
            myc <- constrain(expr)[[pp]]
            if (!is.null(myc)) {
                parval <- mom$parval
                vals <- c(parval,constrainpar,mom$v,mom$e)[attributes(myc)$args]
                fval <- try(myc(unlist(vals)),silent=TRUE)
                fmat <- inherits(fval,"try-error")
                if (fmat) fval <- myc(rbind(unlist(vals)))
                if (!is.null(attributes(fval)$grad)) {
                    if (fmat) {
                        Gr <- attributes(fval)$grad(rbind(unlist(vals)))
                    } else {
                        Gr <- attributes(fval)$grad(unlist(vals))
                    }
                } else {
                    if (fmat) {
                        Gr <- as.numeric(numDeriv::jacobian(myc, rbind(unlist(vals))))
                    } else {
                        Gr <- as.numeric(numDeriv::jacobian(myc, unlist(vals)))
                    }
                }
                mat.idx <- mom$constrain.idx[[pp]]
                cname <- c(cname,pp)
                attributes(fval)$grad <- Gr
                attributes(fval)$vals <- vals
                constrainpar <- c(constrainpar,list(fval)); names(constrainpar) <- cname

                for (jj in seq_len(length(vals))) {
                    allpars <- c(nn$A[attributes(vals[[jj]])$reg.idx[1]],
                                 nn$P[attributes(vals[[jj]])$cov.idx[1]],
                                 nn$v[attributes(vals[[jj]])$m.idx[1]],
                                 nn$e[attributes(vals[[jj]])$e.idx[1]]
                                 )
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
        ii$Ik <- diag(nrow=nobs)
        ii$Im <- diag(nrow=ncol(ii$A))
        ##    ii$Kkk <- commutation(nobs,sparse=FALSE)
    }

    K <- nobs
    ## if (N>10) {
    if (!lava.options()$devel) {
        dG <- with(mom, kronprod(t(IAi),G,res$dA))
        G3 <- with(mom, kronprod(G,G,res$dP))
        GP <- with(mom,G%*%P)
        G1 <- with(mom, kronprod(GP,ii$Ik,dG))
        G2 <- G1[as.vector(matrix(seq_len(K^2),K,byrow=TRUE)),]
        dS <- G1+G2+G3
    } else {
        dG <- with(mom, kronprod(t(IAi),G,res$dA[,ii$parBelongsTo$reg,drop=FALSE]))
        G3 <- with(mom, kronprod(G,G,res$dP[,ii$parBelongsTo$cov,drop=FALSE]))
        GP <- with(mom,G%*%P)
        G1 <- with(mom, kronprod(GP,ii$Ik,dG))
        G2 <- G1[as.vector(matrix(seq_len(K^2),K,byrow=TRUE)),]
        dS <- matrix(0,nrow=nrow(G1),ncol=ncol(res$dA))
        dS[,ii$parBelongsTo$reg] <- G1+G2;  dS[,ii$parBelongsTo$cov] <- G3
    }

    res <- c(res, list(dG=dG, dS=dS))

    if (!is.null(mom$v)) {
        if (lava.options()$devel) {
            dG <- with(mom, kronprod(t(IAi),G,res$dA[,with(ii$parBelongsTo,c(mean,reg)),drop=FALSE]))
        }
        dxi <-
            with(mom, kronprod(rbind(v),dG))
        ##  with(mom, kronprod(rbind(v),ii$Ik,dG))
        if (is.matrix(mom$v) && nrow(mom$v)>1) {
            ## reorder
            k <- nrow(dxi)/nrow(mom$v)
            idx0 <- seq(nrow(mom$v))*k-k+1
            idx <- unlist(lapply(1:k,function(x) idx0+x-1))
            dxi <- dxi[idx,,drop=FALSE]
        }

        if (!is.null(res$dv)) {
            if (!(lava.options()$devel)) {
                if (is.matrix(mom$v) && nrow(mom$v)>1) {
                    dxi <- dxi + (mom$G%*%res$dv)%x%cbind(rep(1,nrow(mom$v)))
                } else {
                    dxi <- dxi+ mom$G%*%res$dv
                }
            } else {
                dxi <- dxi+ mom$G%*%res$dv[,with(ii$parBelongsTo,c(mean,reg))]
            }
        }
        res <- c(res, list(dxi=dxi))
        if (!is.null(mu)) {
            muv <- rbind(mu-mom$xi)
            dT <- suppressMessages(-t(ii$Ik%x%muv + muv%x%ii$Ik) %*% dxi)
            res <- c(res, list(dT=dT))
        }
    }


    if (second) {
        k <- nrow(ii$A)
        K <- ii$Kkk ## commutation(k,k)
        I <- ii$Ik ## diag(k)
        I2 <- diag(nrow=k*k)
        d2S1 <-  t(
        (I %x% K %x% I) %*% (
            ( I2 %x% as.vector(mom$G) )%*% dG +
            ( as.vector(mom$P) %x% I2 )%*% (dP)
        ) %*% t(dG)
        )
        d2S2 <- K%*%d2S1
        d2S3 <- t(
        (I %x% K %x% I) %*% (
            ( I2 %x% as.vector(mom$G) )%*% dG +
            ( as.vector(mom$G) %x% I2 )%*% dG
        ) %*% t(dP)
        )
        vec.d2S <- d2S1+d2S2+d2S3
        res <- c(res, list(d2vecS=vec.d2S))
    }

    return(res)
}
