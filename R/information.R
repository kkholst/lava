##' @export
`information` <-
    function(x,...) UseMethod("information")

###{{{ information.lvm

##' @export
information.lvm <- function(x,p,n,type=ifelse(model=="gaussian",
                                       c("E","hessian","varS","outer","sandwich","robust","num"),"outer"),
                     data,weights=NULL,
                     data2=NULL,
                     model="gaussian",
                     method=lava.options()$Dmethod,
                     inverse=FALSE, pinv=TRUE,
                     score=TRUE,...) {
    if (missing(n))
        n <- NROW(data)
    if (type[1]%in%c("sandwich","robust")) {
        cl <- match.call()
        cl$inverse <- !inverse
        cl$type <- "outer"
        A <- eval.parent(cl)
        cl$inverse <- !(cl$inverse)
        cl$type <- ifelse(type[1]=="sandwich","E","hessian")
        B <- eval.parent(cl)
        return(B%*%A%*%B)
    }
    if (type[1]%in%c("num","hessian","obs","observed")  | (type[1]%in%c("E","hessian") & model!="gaussian")) {
        myf <- function(p0) score(x, p=p0, model=model,data=data, weights=weights,data2=data2,indiv=FALSE,n=n) ##...)
        I <- -numDeriv::jacobian(myf,p,method=method)
        res <- (I+t(I))/2 # Symmetric result
        if (inverse) {
            if (pinv)
                iI <- Inverse(res)
            else
                iI <- solve(res)
            return(iI)
        }
        return(res)
    }
    if (type[1]=="varS" | type[1]=="outer") {
        S <- score(x,p=p,data=na.omit(data),model=model,weights=weights,data2=data2,indiv=TRUE,...)
        res <- t(S)%*%S
        if (inverse) {
            if (pinv)
                iI <- Inverse(res)
            else
                iI <- solve(res)
            return(iI)
        }
        attributes(res)$grad <- colSums(S)
        return(res)
    }

    if (n>1) {
        xfix <- colnames(data)[(colnames(data)%in%parlabels(x,exo=TRUE))]
        xconstrain <- intersect(unlist(lapply(constrain(x),function(z) attributes(z)$args)),manifest(x))

        if (length(xfix)>0 | length(xconstrain)>0) { ##### Random slopes!
            x0 <- x
            if (length(xfix)>0) {
                nrow <- length(vars(x))
                xpos <- lapply(xfix,function(y) which(regfix(x)$labels==y))
                colpos <- lapply(xpos, function(y) ceiling(y/nrow))
                rowpos <- lapply(xpos, function(y) (y-1)%%nrow+1)
                myfix <- list(var=xfix, col=colpos, row=rowpos)
                for (i in seq_along(myfix$var))
                    for (j in seq_along(myfix$col[[i]]))
                        regfix(x0, from=vars(x0)[myfix$row[[i]]][j],to=vars(x0)[myfix$col[[i]]][j]) <-
                            data[1,myfix$var[[i]]]
                index(x0) <- reindex(x0,zeroones=TRUE,deriv=TRUE)
            }
            pp <- modelPar(x0,p)
            p0 <- with(pp, c(meanpar,p,p2))
            myfun <- function(ii) {
                if (length(xfix)>0)
                    for (i in seq_along(myfix$var)) {
                        for (j in seq_along(myfix$col[[i]])) {
                            index(x0)$A[cbind(myfix$row[[i]],myfix$col[[i]])] <- data[ii,myfix$var[[i]]]
                        }
                    }
                ww <- NULL
                if (!is.null(weights))
                    ww <- weights[ii,]
                return(information(x0,p=p,n=1,type=type,weights=ww,data=data[ii,]))
            }
            L <- lapply(seq_len(nrow(data)),function(y) myfun(y))
            val <- apply(array(unlist(L),dim=c(length(p0),length(p0),nrow(data))),c(1,2),sum)
            if (inverse) {
                if (pinv)
                    iI <- Inverse(val)
                else
                    iI <- solve(val)
                return(iI)
            }
            return(val)
        }
    }

    if (!is.null(weights) && is.matrix(weights)) {
        L <- lapply(seq_len(nrow(weights)),function(y) information(x,p=p,n=1,type=type,weights=weights[y,]))
        val <- apply(array(unlist(L),dim=c(length(p),length(p),nrow(weights))),c(1,2),sum)
        if (inverse) {
            if (pinv)
                iI <- Inverse(val)
            else
                iI <- solve(val)
            return(iI)
        }
        return(val)
    }
    mp <- moments(x,p,data=data)
    pp <- modelPar(x,p)
    D <- deriv.lvm(x, meanpar=pp$meanpar, mom=mp, p=p)##, all=length(constrain(x))>0)
    C <- mp$C
    iC <- Inverse(C,det=FALSE, symmetric = TRUE)

    if (is.null(weights)) {
        ##    W <- diag(ncol(iC))
    } else {
        if (length(weights)<ncol(iC)) {
            oldweights <- weights
            weights <- rbind(rep(1,ncol(iC))) ## Ones at exogenous var.
            idx <- index(x)$vars%in%index(x)$exogenous
            print(idx); print(oldweights)
            weights[,idx] <- oldweights
        }
        W <- diag(nrow=as.numeric(weights))
        iW <- W
        diag(iW) <- 1/diag(iW)
    }

    if (is.null(weights)) {
        ## information_Sigma <-  n/2*t(D$dS)%*%((iC)%x%(iC))%*%(D$dS)
        if (lava.options()$devel) {
            information_Sigma <- matrix(0,length(p),length(p))
            imean <- with(index(x)$parBelongsTo,mean)
            information_Sigma[-imean,-imean] <- n/2*t(D$dS[,-imean])%*%kronprod(iC,iC,D$dS[,-imean])
        } else {
            information_Sigma <- n/2*t(D$dS)%*%kronprod(iC,iC,D$dS)
        }
    } else {
        ## information_Sigma <-  n/2*t(D$dS)%*%((iC)%x%(iC%*%W))%*%(D$dS)
        information_Sigma <- n/2*t(D$dS)%*%kronprod(iC,iC%*%W,D$dS)
    }

    dxi <- D$dxi;

    if (!any(dxi>0)) { 
        if (inverse) {
            if (pinv)
                iI <- Inverse(information_Sigma)
            else
                iI <- solve(information_Sigma)
            return(iI)
        }
        return(information_Sigma)
    }
    ii <- index(x)
    if (is.null(weights)) {
        information_mu <- n*t(dxi) %*% (iC) %*% (dxi)
    } else {
        information_mu <- n*t(dxi) %*% (iC%*%W) %*% (dxi)
    }

    if (!(lava.options()$devel)) {
        information <- information_Sigma+information_mu
    } else {
        mparidx <- with(ii$parBelongsTo,c(mean,reg))
        information <- information_Sigma
        information[mparidx,mparidx] <- information[mparidx,mparidx] + information_mu
    }

    if (inverse) {
        if (pinv)
            iI <- Inverse(information, symmetric = TRUE)
        else
            iI <- solve(information)
        return(iI)
    }
    return(information)
}

###}}} information.lvm

###{{{ information.lvmfit

##' @export
information.lvmfit <- function(x,p=pars(x),n=x$data$n,data=model.frame(x),model=x$estimator,weights=Weights(x),
                        data2=x$data$data2,
                        ...) {
    I <- information(x$model0,p=p,n=n,data=data,model=model,
                     weights=weights,data2=data2,...)
    if (ncol(I)<length(p)) {
        I <- blockdiag(I,matrix(0,length(p)-ncol(I),length(p)-ncol(I)))
    }
    return(I)
}

###}}} information.lvmfit


##' @export
information.lvm.missing <- function(x,
                             p=coef(x), estimator=x$estimator,
                             weights=Weights(x$estimate),
                             ...) {
    information(x$estimate$model0, p=p, model=estimator, weights=weights,...)
}

##' @export
information.multigroupfit <- function(x,p=pars(x), weights=Weights(x), estimator=x$estimator, ...) {
    information(x$model0,p=p, weights=weights, model=estimator ,...)
}

##' @export
information.multigroup <- function(x,data=x$data,weights=NULL,p,indiv=FALSE,...) {
    rm <- procrandomslope(x)
    pp <- with(rm, modelPar(model,p)$p)
    parord <- modelPar(rm$model,seq_len(with(rm$model,npar+npar.mean)))$p
    I <- matrix(0,nrow=length(p),ncol=length(p))
    if (!indiv) {
        for (i in seq_len(x$ngroup))
            I[parord[[i]],parord[[i]]] <- I[parord[[i]],parord[[i]]] + information(x$lvm[[i]],p=pp[[i]],data=data[[i]],weights=weights[[i]],...)
    } else {
        I <- list()
        for (i in seq_len(x$ngroup))
            I <- c(I, list(information(x$lvm[[i]],p=pp[[i]],data=data[[i]],weights=weights[[i]],...)))
    }
    return(I)
}
