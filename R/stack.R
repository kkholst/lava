##' Stack estimating equations
##'
##' Stack estimating equations (two-stage estimator)
##' @param x Model 1
##' @param model2 Model 2
##' @param D1u Derivative of score of model 2 w.r.t. parameter vector of model 1
##' @param inv.D2u Inverse of deri
##' @param propensity propensity score (vector or function)
##' @param dpropensity derivative of propensity score wrt parameters of model 1
##' @param U Optional score function (model 2) as function of all parameters 
##' @param keep1 If FALSE only parameters of model 2 is returned
##' @param propensity.arg Arguments to propensity function
##' @param estimate.arg Arguments to 'estimate'
##' @param na.action Method for dealing with missing data in propensity score
##' @param ... Additional arguments to lower level functions
##' @aliases stack.estimate
##' @export
##' @examples
##' m <- lvm(z0~x)
##' Missing(m, z ~ z0) <- r~x
##' distribution(m,~x) <- binomial.lvm()
##' p <- c(r=-1,'r~x'=0.5,'z0~x'=2)
##' beta <- p[3]/2
##' d <- sim(m,500,p=p,seed=1)
##' m1 <- estimate(r~x,data=d,family=binomial)
##' d$w <- d$r/predict(m1,type="response")
##' m2 <- estimate(z~1, weights=w, data=d)
##' (e <- stack(m1,m2,propensity=TRUE))
stack.estimate <- function(x,model2,D1u,inv.D2u,
                   propensity,dpropensity,U,
                   keep1=FALSE,
                   propensity.arg,estimate.arg,
                   na.action=na.pass, ...) {
    ic2 <- IC(model2)
    if (!missing(U)) {
        D1u <- numDeriv::jacobian(U, coef(x))
        tryCatch(D2u <- numDeriv::jacobian(function(p)
          U(coef(x),p), coef(model2)),
          error=function(...) NULL)
        if (!is.null(D2u))
          inv.D2u <- solve(D2u)
    }
    if (missing(inv.D2u)) {
        inv.D2u <- -attributes(ic2)$bread
    }
    if (is.null(inv.D2u)) stop("Need derivative of second stage score")
    if (!missing(propensity) && is.logical(propensity) && propensity) {
        if (!inherits(x,"glm")) stop("'x' must be of 'glm' class")
        f <- family(x)
        ginv <- f$linkinv
        dginv <- f$mu.eta ## D[linkinv]        
        dpropensity <- function(p, data=model.frame(x), ...) {
            X <- model.matrix(formula(x),data, ...)
            p1 <- dginv(X%*%p)
            apply(X,2,function(x) x*p1)
        }
        propensity <- function(p, data=model.frame(x), ...) {
            X <- model.matrix(formula(x),data, ...)
            ginv(X%*%p)            
        }        
    }    
    if (!missing(propensity) && is.function(propensity)) {
        op <- options(na.action=na.action)
        wfun <- propensity
        if (!missing(propensity.arg)) {
            wfun <- function(p,...)
                do.call(propensity, c(list(p),propensity.arg))
        }
        if (missing(dpropensity)) {
            dpropensity <- numDeriv::jacobian(wfun,coef(x),...)
        } else if (is.function(dpropensity)) {
            if (!missing(propensity.arg)) {
                dwfun <- function(p,...)
                    do.call(dpropensity, c(list(p),propensity.arg))
                dpropensity <- dwfun(coef(x),...)
            } else {
                dpropensity <- dpropensity(coef(x),...)
            }
        }
        propensity <- wfun(coef(x),...)
        mis <- which(is.na(propensity))
        propensity[mis] <- 0
        dpropensity[mis,] <- 0
        options(op)
    }
    if (!missing(dpropensity)) {
        D2u <- Inverse(inv.D2u)
        u2 <- -ic2%*%D2u ## Score of stage two equation derived from estimated influence function
        ## Derivative of score wrt first set of parameters (propensity score model)
        D1u <- crossprod(apply(u2,2,function(x) -na.pass0(x/propensity)),dpropensity) / NROW(ic2)
    }
    ii <- IC(merge(x,model2))
    ic1. <- ii[,seq_along(coef(x)),drop=FALSE]
    ic2. <- ii[,length(coef(x))+seq_along(coef(model2)),drop=FALSE]
    ic3 <- -t(inv.D2u%*%(D1u%*%t(ic1.)))
    val <- cbind(ic2. + ic3)
    if (!keep1) return(estimate(coef=coef(model2), IC=val, ...))
    estimate(coef=c(coef(x), coef(model2)),
             IC=cbind(ic1., val)*NROW(val), ...)
}

##' @export
stack.glm <- function(x, model2, ...) {
    stack.estimate(x, model2, ...)
}

##' @export
stack.lvmfit <- function(x, model2, ...) {
    stack(estimate(x), model2, ...)
}

