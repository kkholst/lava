##' @export
stack.estimate <- function(x,y,Ix,iIy,weight,dweight,Ualpha,k=1,...) {
    iid1 <- iid(x)
    iid2 <- iid(y)
    if (missing(iIy)) {
        iIy <- attributes(iid2)$bread
    }
    if (is.null(iIy)) stop("Need derivative of second stage score")
    if (!missing(Ualpha)) {
        Ix <- numDeriv::jacobian(Ualpha,coef(x))
    }
    if (!missing(weight) && is.function(weight)) {
        dweight <- numDeriv::jacobian(weight,coef(x))
        weight <- weight(coef(x))
    }
    if (!missing(dweight)) {
        Iy <- Inverse(iIy,tol=1e-9)
        u2 <- iid2%*%Iy ## Score of stage two equation derived from estimated influence function
        ## Derivative of score wrt first set of parameters (weight-parameters)
        Ix <- crossprod(apply(u2,2,function(x) -x/weight),dweight)
    }
    ii <- iid(merge(x,y))
    iid1. <- ii[,seq_along(coef(x)),drop=FALSE]
    iid2. <- ii[,length(coef(x))+seq_along(coef(y)),drop=FALSE] 
    iid3 <- t(iIy%*%(Ix%*%t(iid1.)))
    estimate(coef=c(coef(x),coef(y)),iid=cbind(iid1.,iid2. + k*iid3))
}

##' @export
stack.glm <- function(x,y,...) {
    stack(estimate(x),estimate(y),...)
}
