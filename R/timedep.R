##' Add time-varying covariate effects to model
##'
##' @title Time-dependent parameters
##' @param object Model
##' @param formula Formula with rhs specifying time-varying covariates
##' @param rate Optional rate parameters. If given as a vector this
##' parameter is interpreted as the raw (baseline-)rates within each
##' time interval defined by \code{timecut}.  If given as a matrix the
##' parameters are interpreted as log-rates (and log-rate-ratios for
##' the time-varying covariates defined in the formula).
##' @param timecut Time intervals
##' @param type Type of model (default piecewise constant intensity)
##' @param ... Additional arguments to lower level functions
##' @author Klaus K. Holst
##' @aliases timedep timedep<-
##' @export
##' @examples
##'
##' ## Piecewise constant hazard
##' m <- lvm(y~1)
##' m <- timedep(m,y~1,timecut=c(0,5),rate=c(0.5,0.3))
##'
##' \dontrun{
##' d <- sim(m,1e4); d$status <- TRUE
##' dd <- mets::lifetable(Surv(y,status)~1,data=d,breaks=c(0,5,10));
##' exp(coef(glm(events ~ offset(log(atrisk)) + -1 + interval, dd, family=poisson)))
##' }
##'
##'
##' ## Piecewise constant hazard and time-varying effect of z1
##' m <- lvm(y~1)
##' distribution(m,~z1) <- ones.lvm(0.5)
##' R <- log(cbind(c(0.2,0.7,0.9),c(0.5,0.3,0.3)))
##' m <- timedep(m,y~z1,timecut=c(0,3,5),rate=R)
##'
##' \dontrun{
##' d <- sim(m,1e4); d$status <- TRUE
##' dd <- mets::lifetable(Surv(y,status)~z1,data=d,breaks=c(0,3,5,Inf));
##' exp(coef(glm(events ~ offset(log(atrisk)) + -1 + interval+z1:interval, dd, family=poisson)))
##' }
##'
##'
##'
##' ## Explicit simulation of time-varying effects
##' m <- lvm(y~1)
##' distribution(m,~z1) <- ones.lvm(0.5)
##' distribution(m,~z2) <- binomial.lvm(p=0.5)
##' #variance(m,~m1+m2) <- 0
##' #regression(m,m1[m1:0] ~ z1) <- log(0.5)
##' #regression(m,m2[m2:0] ~ z1) <- log(0.3)
##' regression(m,m1 ~ z1,variance=0) <- log(0.5)
##' regression(m,m2 ~ z1,variance=0) <- log(0.3)
##' intercept(m,~m1+m2) <- c(-0.5,0)
##' m <- timedep(m,y~m1+m2,timecut=c(0,5))
##'
##' \dontrun{
##' d <- sim(m,1e5); d$status <- TRUE
##' dd <- mets::lifetable(Surv(y,status)~z1,data=d,breaks=c(0,5,Inf))
##' exp(coef(glm(events ~ offset(log(atrisk)) + -1 + interval + interval:z1, dd, family=poisson)))
##' }
timedep <- function(object,formula,rate,timecut,type="coxExponential.lvm",...) {
    if (missing(timecut)) stop("'timecut' needed")
    ##if (inherits(formula,"formula"))
    ff <- getoutcome(formula)
    simvars <- attributes(ff)$x
    if (is.null(object$attributes$simvar)) {
        object$attributes$simvar <- list(simvars)
        names(object$attributes$simvar) <- ff
        object$attributes$timedep <- object$attributes$simvar
    } else {
        object$attributes$simvar[[ff]] <- simvars
        object$attributes$timedep[[ff]] <- simvars
    }
    if (missing(rate)) rate <- rep(1,length(timecut))
    
    args <- list(timecut=timecut,rate=rate,...)
    covariance(object,ff) <- 1
    distribution(object,ff) <- do.call(type,args)
    return(object)
}

##' @export
"timedep<-" <- function(object,...,value) {
    timedep(object,value,...)
}
