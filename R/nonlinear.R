##' @export
"nonlinear<-" <- function(object,...,value) UseMethod("nonlinear<-")

##' @export
"nonlinear" <- function(object,...) UseMethod("nonlinear")

naturalcubicspline <- function(x, knots=stats::median(x,na.rm=TRUE), boundary=range(x,na.rm=TRUE)) {
    ## C2 functions, piecewise cubic
    breaks <- c(boundary[1],knots,boundary[2])
    K <- length(breaks)
    g <- function(x,tau) (x-tau)^3*((x-tau)>0)
    gg <- matrix(0,nrow=length(x),ncol=K)
    for (i in seq(K)) {
        gg[,i] <- g(x,breaks[i])
    }
    B <- matrix(0,nrow=length(x),ncol=K-2)
    for (i in seq(K-2)) {
        B[,i] <- gg[,i] -
            (breaks[K]-breaks[i])/(breaks[K]-breaks[K-1])*gg[,K-1] +
            (breaks[K-1]-breaks[i])/(breaks[K]-breaks[K-1])*gg[,K]
    }
    cbind(x,B)
}

ncspred <- function(mu, var, knots=c(-5,0,5)) {
    breaks <- knots
    K <- length(breaks)

    v <- as.vector(var)
    k <- sqrt(v/(2*pi))
    g <- function(x,tau) {
        x0 <- (x-tau)
        x2 <- x0^2
        p0 <- 1-pnorm(-x0/sqrt(v)) # P(x>tau|...)
        k*(2*v + x2)*exp(-(x0/(sqrt(2*v)))^2) +
            x0*(x2+3*v)*p0
    }
    n <- NROW(mu)
    gg <- matrix(0,nrow=n,ncol=K)
    for (i in seq(K)) {
        gg[,i] <- g(mu,breaks[i])
    }
    B <- matrix(0,nrow=n,ncol=K-2)
    for (i in seq(K-2)) {
        B[,i] <- gg[,i] -
            (breaks[K]-breaks[i])/(breaks[K]-breaks[K-1])*gg[,K-1] +
            (breaks[K-1]-breaks[i])/(breaks[K]-breaks[K-1])*gg[,K]
    }
    cbind(mu,B)
}


##' @export
nonlinear.lvm <- function(object, to, from=NULL, type=c("quadratic"), knots=c(-5,0,5), names, ...) {
    if (missing(to)) {
        return(object$attributes$nonlinear)
    }
    if (inherits(to,"formula")) {
        yy <- decomp.specials(getoutcome(to))
        myvars <- all.vars(to)
        from <- setdiff(myvars,yy)
        to <- yy
    }
    if (length(to)>1) stop("Supply only one response variable")
    if (length(from)>1) stop("Supply only one explanatory variable")
    object <- cancel(object, c(from,to))
    variance(object) <- to
    f <- pred <- NULL

    if (tolower(type)[1]%in%c("ncs","spline","naturalspline","cubicspline","natural cubic spline")) {
        if (is.null(knots)) stop("Need cut-points ('knots')")
        if (length(knots)<3) {
            warning("Supply at least three knots (one interior and boundaries)")
            ## Fall-back to linear
            type <- "linear"
        }
        if (missing(names)) names <- paste0(from,"_",seq(length(knots)-1))
        f <- function(p,x) {
            B <- cbind(1,naturalcubicspline(x,knots=knots[-c(1,length(knots))],boundary=knots[c(1,length(knots))]))
            colnames(B) <- c("(Intercept)",names)
            as.vector(B%*%p)
        }
        pred <- function(mu,var,...) {
            B <- ncspred(mu,var,knots=knots)
            structure(B,dimnames=list(NULL,names))
        }
    }

    
    if (tolower(type)[1]=="linear") {
        if (missing(names)) names <- from
        f <- function(p,x) p[1] + p[2]*x
        pred <- function(mu,var,...) {            
            structure(cbind(mu[,1]),dimnames=list(NULL,names))
        }
    }
    
    if (tolower(type)[1]=="quadratic") {
        if (missing(names)) names <- paste0(from,"_",1:2)
        f <- function(p,x) p[1] + p[2]*x + p[3]*(x*x)
        pred <- function(mu,var,...) {
            structure(cbind(mu[,1],mu[,1]^2+var[1]),dimnames=list(NULL,names))
        }
    }

    if (tolower(type)[1]%in%c("piecewise","piecewise linear","linear")) {
        if (is.null(knots)) stop("Need cut-points ('knots')")
    }

    if (tolower(type)[1]%in%c("exp","exponential")) {
        if (missing(names)) names <- paste0(from,"_",1:2)
        f <- function(p,x) p[1] + p[2]*x + p[3]*exp(x)
        pred <- function(mu,var,...) {
            structure(cbind(mu[,1], exp(0.5*var[1] + mu[,1])),dimnames=list(NULL,names))
        }
    }

    if (tolower(type)[1]%in%c("exp0")) {
        if (missing(names)) names <- paste0(from,"_",1)
        f <- function(p,x) p[1] + p[2]*exp(x)
        pred <- function(mu,var,...) {
            structure(cbind(exp(0.5*var[1] + mu[,1])),dimnames=list(NULL,names))
        }
    }

    object$attributes$nonlinear[[to]] <- list(x=from, p=length(names)+1, newx=names, f=f, pred=pred, type=tolower(type[1]))
    return(object)
}

##' @export
nonlinear.lvmfit <- function(object, to, ...) {
    if (missing(to)) {
        return(Model(object)$attributes$nonlinear)
    }
    Model(object) <- nonlinear(Model(object),to=to,...)
    return(object)
}

##' @export
nonlinear.twostage.lvm <- function(object, ...) {
    return(object$nonlinear)
}

##' @export
nonlinear.lvmfit <- function(object, ...) {
    return(object$nonlinear)
}

##' @export
`nonlinear<-.lvm` <- function(object, ..., type="quadratic", value) {
    nonlinear(object,to=value,type=type,...)
}
