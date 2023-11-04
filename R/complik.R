##' Composite Likelihood for probit latent variable models
##'
##' Estimate parameters in a probit latent variable model via a composite
##' likelihood decomposition.
##' @param x \code{lvm}-object
##' @param data data.frame
##' @param k Size of composite groups
##' @param type Determines number of groups. With \code{type="nearest"} (default)
##' only neighboring items will be grouped, e.g. for \code{k=2}
##' (y1,y2),(y2,y3),... With \code{type="all"} all combinations of size \code{k}
##' are included
##' @param pairlist A list of indices specifying the composite groups. Optional
##' argument which overrides \code{k} and \code{type} but gives complete
##' flexibility in the specification of the composite likelihood
##' @param messages Control amount of messages printed
##' @param estimator Model (pseudo-likelihood) to use for the pairs/groups
##' @param quick If TRUE the parameter estimates are calculated but all additional
##' information such as standard errors are skipped
##' @param \dots Additional arguments parsed on to lower-level functions
##' @return An object of class \code{estimate.complik} inheriting methods from \code{lvm}
##' @author Klaus K. Holst
##' @seealso estimate
##' @keywords models regression
##' @export
##' @examples
##' m <- lvm(c(y1,y2,y3)~b*x+1*u[0],latent=~u)
##' ordinal(m,K=2) <- ~y1+y2+y3
##' d <- sim(m,50,seed=1)
##' if (requireNamespace("mets", quietly=TRUE)) {
##'    e1 <- complik(m,d,control=list(trace=1),type="all")
##' }
complik <- function(x, data, k=2, type=c("all","nearest"), pairlist,
             messages=0, estimator="normal", quick=FALSE, ...) {

    y <- setdiff(endogenous(x),latent(x))
    binsurv <- rep(FALSE,length(y))
    for (i in 1:length(y)) {
        z <- try(data[,y[i]],silent=TRUE)
        ## binsurv[i] <- is.Surv(z) | (is.factor(z) && length(levels(z))==2)
        if (!inherits(z,"try-error"))
            binsurv[i] <- inherits(z,"Surv") | (is.factor(z))
    }

    ord <- ordinal(x)
    binsurv <- unique(c(y[binsurv],ord))
    if (!missing(pairlist)) {
        binsurvpos <- which(colnames(data)%in%endogenous(x))
    } else {
        binsurvpos <- which(colnames(data)%in%binsurv)
    }
    if (missing(pairlist)) {
        if (type[1]=="all") {
            mypar <- combn(length(binsurv),k) ## all pairs (or multiplets), k=2: k*(k-1)/2
        } else {
            mypar <- sapply(0:(length(binsurv)-k), function(x) x+1:k)
        }
    } else {
        mypar <- pairlist
    }

    if (is.matrix(mypar)) {
        mypar0 <- mypar; mypar <- c()
        for (i in seq(ncol(mypar0)))
            mypar <- c(mypar, list(mypar0[,i]))
    }
    nblocks <- length(mypar)

    mydata0 <- data[,,drop=FALSE]
    mydata <-  as.data.frame(matrix(NA, nblocks*nrow(data), ncol=ncol(data)))
    names(mydata) <- names(mydata0)
    for (i in 1:ncol(mydata)) {
        if (is.factor(data[,i])) {
            mydata[,i] <- factor(mydata[,i],levels=levels(mydata0[,i]))
        }
        if (survival::is.Surv(data[,i])) {
            S <- data[,i]
            for (j in 2:nblocks) S <- rbind(S,data[,i])
            S[] <- NA
            mydata[,i] <- S
        }
    }
    for (ii in 1:nblocks) {
        data0 <- data;
        for (i in binsurvpos[-mypar[[ii]]]) {
            if (survival::is.Surv(data[,i])) {
                S <- data0[,i]
                S[] <- NA
                data0[,i] <- S
            } else {
                data0[,i] <- NA
                if (is.factor(data[,i])) data0[,i] <- factor(data0[,i],levels=levels(data[,i]))
            }
        }
        mydata[(1:nrow(data))+(ii-1)*nrow(data),] <- data0
    }
    suppressWarnings(e0 <- estimate(x,data=mydata,estimator=estimator,missing=TRUE,messages=messages,
                                    hessian=!quick,
                                    ...))
    e0$n <- nrow(data)
    e0$eigenvalues <- rep(1,nrow(e0$coef))
    e0$coef[,2:4] <- NA
    if (!quick) {
        S <- score(e0,indiv=TRUE)
        nd <- nrow(data)
        block1 <- which((1:nd)%in%(rownames(S)))
        blocks <- sapply(1:nblocks, function(x) 1:length(block1)+length(block1)*(x-1))
        if (nblocks==1) {
            Siid <- S
        } else {
            Siid <- matrix(0,nrow=length(block1),ncol=ncol(S))
            for (j in 1:ncol(blocks)) {
                Siid <- Siid+S[blocks[,j],]
            }
        }
        ##B <- solve(information(e0, type="hessian"))
        D <- numDeriv::jacobian(function(p) score(e0, p=p), coef(e0), method=lava.options()$Dmethod)
        B <- -Inverse(D)*NROW(Siid)
        A <- crossprod(Siid)/NROW(Siid)
        e0$bread <- B
        e0$meat <- A
        e0$iidscore <- Siid
        e0$blocks <- blocks
        BA <- B%*%A
        eig <- eigen(BA)$values
        e0$eigenvalues <- eig
        e0$vcov <- (BA%*%B)/NROW(Siid)
        cc <- e0$coef;
        cc[,2] <- sqrt(diag(e0$vcov))
        cc[,3] <- cc[,1]/cc[,2]; cc[,4] <- 2*(1-pnorm(abs(cc[,3])))
        e0$coef <- cc
    }
    class(e0) <- c("estimate.complik",class(e0))
    return(e0)
}

##' @export
score.estimate.complik <- function(x,indiv=FALSE,...) {
    if (!indiv)
        return(colSums(x$iidscore))
    x$iidscore
}

##' @export
IC.estimate.complik <- function(x,...) {
    IC.default(x,bread=x$bread,...)
}

##' @export
logLik.estimate.complik <- function(object,...) {
    with(object, structure(opt$objective,
              nall=n,
              nobs=n,
              df=sum(eigenvalues),
              class="logLik"))
}
