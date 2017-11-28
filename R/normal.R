introotpn <- function(p) {
    ## Find integer root of x^2-x-2*p=0
    n <- 0.5*(1+sqrt(1+8*p))
    if (floor(n)!=n) n <- NA
    return(n)
}
rho2sigma <- function(rho) {
    if (length(rho)==1) return(diag(2)*(1-rho)+rho)
    p <- introotpn(length(rho))
    if (is.na(p)) stop("Unexpected length of correlation coefficients (p=n*(n-1)/2).")
    sigma <- diag(nrow=p)
    offdiag(sigma,type=2) <- rho
    offdiag(sigma,type=3) <- offdiag(t(sigma),type=3)
    return(sigma)
}

##' @export
rmvn <- function(n,mu,sigma,rho,...) {
    if (!missing(rho)) sigma <- rho2sigma(rho)
    if (!missing(mu) && missing(sigma)) sigma <- diag(nrow=length(mu))
    if (missing(sigma)) sigma <- matrix(1)
    if (is.vector(sigma)) sigma <- diag(sigma,ncol=length(sigma))
    if (missing(mu)) mu <- rep(0,ncol(sigma))    
    PP <- with(svd(sigma), v%*%diag(sqrt(d),ncol=length(d))%*%t(u))
    res <- matrix(rnorm(ncol(sigma)*n),ncol=ncol(sigma))%*%PP
    if (NROW(mu)==nrow(res) && NCOL(mu)==ncol(res)) return(res+mu)
    return(res+cbind(rep(1,n))%*%mu)
}

##' @export
dmvn <- function(x,mu,sigma,rho,log=FALSE,nan.zero=TRUE,norm=TRUE,...) {
    if (!missing(rho)) sigma <- rho2sigma(rho)
    if (!missing(mu) && missing(sigma)) sigma <- diag(nrow=length(mu))
    if (missing(sigma)) sigma <- matrix(1)
    if (is.vector(sigma)) sigma <- diag(sigma,ncol=length(sigma))
    if (missing(mu)) mu <- rep(0,ncol(sigma))
    
    if (length(sigma)==1) {
        k <- 1
        isigma <- structure(cbind(1/sigma),det=as.vector(sigma))

    } else {
        k <- ncol(sigma)
        isigma <- Inverse(sigma)
    }
    if (!missing(mu)) {
        if (NROW(mu)==NROW(x) && NCOL(mu)==NCOL(x)) {
            x <- x-mu
        } else {
            x <- t(t(x)-mu)
        }
    }
    logval <- -0.5*(base::log(2*base::pi)*k+
                    base::log(attributes(isigma)$det)+
                    rowSums((x%*%isigma)*x))
    if (nan.zero) logval[is.nan(logval)] <- -Inf
    if (log) return(logval)
    return(exp(logval))
}


normal_method.lvm <- "nlminb0"

normal_objective.lvm <- function(x,p,data,weights=NULL,data2=NULL,indiv=FALSE,...) {
    if (!requireNamespace("mets",quietly=TRUE)) stop("'mets' package required")
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    save.seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", save.seed, envir = .GlobalEnv))
    set.seed(1)
    ii <- lava::index(x)
    y.idx <- ii$endo.idx
    x.idx <- ii$exo.idx
    y <- ii$endogenous
    ord <- lava::ordinal(x)
    atr <- attributes(ord)
    ord <- intersect(y,ord)
    attributes(ord) <- atr

    status <- rep(0,length(y))
    status[match(ord,y)] <- 2

    Table <- (length(y)==length(ord)) && (length(x.idx)==0)
    if (Table) {
        pat <- mets::fast.pattern(data[,y,drop=FALSE],categories=max(data[,y,drop=FALSE])+1)
        data <- pat$pattern
        colnames(data) <- y
    }

    mu <- predict(x,data=data,p=p)
    S <- attributes(mu)$cond.var
    class(mu) <- "matrix"
    thres <- matrix(0,nrow=length(y),max(1,attributes(ord)$K-1)); rownames(thres) <- y
    for (i in seq_len(length(attributes(ord)$fix))) {
        nn <- names(attributes(ord)$idx)[i]
        ii <- attributes(ord)$idx[[nn]]
        val <- (attributes(mu)$e[ii])
        thres[nn,seq_len(length(val))] <-
            cumsum(c(val[1],exp(val[-1])))
    }

    yl <- yu <- as.matrix(data[,y,drop=FALSE])
    if (!inherits(yl[1,1],c("numeric","integer","logical")) ||
        !inherits(yu[1,1],c("numeric","integer","logical")))
        stop("Unexpected data (normal_objective)")

    if (!is.null(data2)) {
        yu[,colnames(data2)] <- data2
        status[match(colnames(data2),y)] <- 1
    }

    l <- mets::loglikMVN(yl,yu,status,mu,S,thres)
    
    if (!is.null(weights)) {
        ##if (is.matrix(weights)) weights <- weights[,1]
        l <- l*weights
    }

    if (Table) {
        l <- l[pat$group+1]
    }
    if (indiv) return(-l)
    return(-sum(l))
}


normal_logLik.lvm <- function(object,p,data,data2=NULL,...) {
    res <- -normal_objective.lvm(x=object,p=p,data=data,data2=data2,...)
    return(res)
}

normal_gradient.lvm <- function(x,p,data,weights=NULL,data2=NULL,indiv=FALSE,...) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    save.seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", save.seed, envir = .GlobalEnv))
    if (!requireNamespace("mets",quietly=TRUE)) stop("'mets' package required")
    if  (is.null(ordinal(x)) && is.null(data2) && is.null(weights)) {
        D <- deriv.lvm(x,p=p)
        M <- moments(x,p)
        Y <- as.matrix(data[,manifest(x)])
        mu <- M$xi%x%rep(1,nrow(Y))
        ss <- -mets::scoreMVN(Y,mu,M$C,D$dxi,D$dS)
        if (!indiv) return(colSums(ss))
        return(ss)
    }
    if (indiv) {
        return(numDeriv::jacobian(function(p0) normal_objective.lvm(x,p=p0,data=data,weights=weights,data2=data2,indiv=TRUE,...),p,method=lava.options()$Dmethod))
    }
    numDeriv::grad(function(p0) normal_objective.lvm(x,p=p0,data=data,weights=weights,data2=data2,...),p,method=lava.options()$Dmethod)
}

normal_hessian.lvm <- function(x,p,outer=FALSE,data2=NULL,...) {
    if (!requireNamespace("mets",quietly=TRUE)) stop("'mets' package required")
    dots <- list(...); dots$weights <- NULL
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    save.seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", save.seed, envir = .GlobalEnv))
    if (!outer) {
        f <- function(p) {
            set.seed(1)
        do.call("normal_objective.lvm", c(list(x,p=p,indiv=FALSE,data2=data2),dots))
        }
        g <- function(p) {
            set.seed(1)
            do.call("normal_gradient.lvm", c(list(x,p=p,indiv=FALSE,data2=data2),dots))
        }
        if (is.null(ordinal(x)) && is.null(data2))
            return(numDeriv::jacobian(g,p))
        else {
            return(numDeriv::hessian(f,p))
        }
    }
    ## Else calculate outer product of the score (empirical variance of score)
    S <- normal_gradient.lvm(x,p=p,indiv=TRUE,...)
    J <- t(S)%*%S
    attributes(J)$grad <- colSums(S)
    return(J)

}

##normal_gradient.lvm <- normal_hessian.lvm <- NULL
