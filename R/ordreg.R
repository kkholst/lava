ordreg_threshold <- function(theta) {
    v <- theta[1]
    if (length(theta)>1) v <- cumsum(c(v,exp(theta[seq(length(theta)-1L)+1L])))
    return(v)
}

ordreg_ithreshold <- function(v) {
    theta <- v[1]
    if (length(v)>1) theta <- c(theta,log(-rev(diff(rev(v)))))
    return(theta)
}

ordreg_dthreshold <- function(theta) {
    K <- length(theta)+1
    Da <- matrix(0,K,K-1)
    Da[seq(K-1),1L] <- 1L
    for (i in seq_len(K-2)+1) Da[seq(i,K-1),i] <- exp(theta[i])
    Da
}

##' Ordinal regression models
##'
##' @title Univariate cumulative link regression models
##' @param formula formula
##' @param data data.frame
##' @param offset offset
##' @param family family (default proportional odds)
##' @param start optional starting values
##' @param fast If TRUE standard errors etc. will not be calculated
##' @param ... Additional arguments to lower level functions
##' @details Let \eqn{Y\in\{1,...,J\}} be the ordinal outcome and \eqn{X} a
##'   vector of covariates. The cumulative link model is given by \deqn{ P(Y\leq
##'   j|X=x) = g(a_j - b^\top x), j=1,...,J-1.} The default link function is the
##'   Probit function, i.e. where \eqn{g} is equal to the standard normal cumulative
##'   distribution function. The proportional odds model is obtained with
##'   \code{family=binomial(logit)}.
##'
##'   Note, the intercept parameters are parametrized such that they are
##'   monotone increasing \eqn{a_1 < \cdots < a_{J-1}}. To get the parameter
##'   estimates of the actual \eqn{a_j}'s use the \code{summary} method.
##' @export
##' @author Klaus K. Holst
##' @examples
##' m <- lvm(y~x)
##' ordinal(m,K=3) <- ~y
##' d <- sim(m,100)
##' e <- ordreg(y~x,d)
##' summary(e)
ordreg <- function(formula,data=parent.frame(),offset,family=stats::binomial("probit"),start,fast=FALSE,...) {
    y <- ordered(model.frame(update(formula,.~0),data)[,1])
    lev <- levels(y)
    X <- model.matrix(update(formula,.~.+1),data=data)[,-1,drop=FALSE]
    up <- new.env()
    assign("h",family$linkinv,envir=up)
    assign("dh",family$mu.eta,envir=up)
    assign("y",as.numeric(y),envir=up)
    assign("X",X,envir=up)
    assign("K",nlevels(y),envir=up)
    assign("n",length(y),envir=up)
    assign("p",NCOL(X),envir=up)
    assign("threshold",
           function(theta,K) ordreg_threshold(theta[seq(K-1)]), envir=up)
    assign("dthreshold",
           function(theta,K) ordreg_dthreshold(theta[seq(K-1)]), envir=up)
    ff <- function(theta) -ordreg_logL(theta,up)
    gg <- function(theta) -ordreg_score(theta,up)
    if (missing(start)) start <- with(up,c(rep(-1,up$K-1),rep(0,p)))
    op <- nlminb(start,ff,gg)
    cc <- op$par;
    if (fast) return(structure(cc,threshold=up$threshold(cc,up$K))) ##,up$K)))
    nn <- c(paste(lev[-length(lev)], lev[-1L], sep = "|"),
                   colnames(X))
    I <- -ordreg_hessian(cc,up)
    names(cc) <- nn
    dimnames(I) <- list(nn,nn)
    res <- list(vcov=solve(I),coef=cc,call=match.call(),up=up,opt=op,levels=lev)
    structure(res,class="ordreg")
}

##' @export
print.ordreg <- function(x,...) {
    cat("Call:\n"); print(x$call)
    cat("\nParameter Estimates:\n")
    print(x$coef)
}

##' @export
summary.ordreg <- function(object, contrast = NULL, ...) {
  lev <- object$levels
  J <- length(lev)
  est <- summary(estimate(object, function(x) c(ordreg_threshold(x[seq_len(J-1)]),
                                                tail(x, length(x)-J+1))),
                 contrast = contrast,...)
  est$call <- NULL
  res <- list(coef=est, logLik=logLik(object), AIC=AIC(object))
  class(res) <- "summary.ordreg"
  return(res)
}

##' @export
print.summary.ordreg <- function(x,alpha=0.95,...) {
  cat("AIC: ", x$AIC, "\n\n")
  print(x$coef)
  cat("\n")
}

##' @export
score.ordreg <- function(x,p=coef(x),indiv=FALSE,...) {
    ordreg_score(p,x$up)
    if (!indiv) return(colSums(x$up$score))
    x$up$score
}

##' @export
logLik.ordreg <- function(object,p=coef(object),indiv=FALSE,...) {
    ordreg_logL(p,object$up)
    res <- log(object$up$pr)
    if (!indiv) res <- sum(res)
    structure(res,nall=length(object$up$pr),nobs=object$up$pr,df=length(p),class="logLik")
}

##' @export
coef.ordreg <- function(object,...) object$coef

##' @export
vcov.ordreg <- function(object,...) object$vcov

ordreg_logL <- function(theta,env,indiv=FALSE,...) {
    if (length(theta)!=with(env,p+K-1)) stop("Wrong dimension")
    env$theta <- theta
    if (env$p>0) beta <- with(env,theta[seq(p)+K-1])
    alpha <- with(env, threshold(theta,K))
    env$alpha <- alpha
    env$beta <- beta
    if (env$p>0) eta <- env$X%*%beta else eta <- cbind(rep(0,env$n))
    env$lp <- kronecker(-eta,rbind(alpha),"+")
    F <- with(env,h(lp))
    Pr <- cbind(F,1)-cbind(0,F)
    pr <- Pr[with(env,cbind(seq(n),as.numeric(y)))]
    env$pr <- pr
    sum(log(pr))
}

ordreg_score <- function(theta,env,...) {
    if (!identical(theta,env$theta)) ordreg_logL(theta,env)
    Da <- with(env,dthreshold(theta,K))
    dF <- with(env, cbind(dh(lp),0))
    idx1 <- with(env,which(as.numeric(y)==1))
    S1 <- cbind(Da[as.numeric(env$y),,drop=FALSE],-env$X)
    S1 <- dF[with(env,cbind(seq(n),as.numeric(y)))]*S1
    y2 <- env$y-1; y2[idx1] <- env$K
    S2 <- cbind(Da[y2,,drop=FALSE],-env$X)
    S2 <- dF[cbind(seq(env$n),y2)]*S2
    env$score <- 1/env$pr*(S1-S2)
    colSums(env$score)
}
ordreg_hessian <- function(theta,env,...) {
    numDeriv::jacobian(function(p) ordreg_score(p,env,...),theta,...)
}

##' @export
predict.ordreg <- function(object,p=coef(object),type=c("prob","cumulative"),...) {
    env <- object$up
    env$theta <- p
    if (env$p>0) beta <- with(env,theta[seq(p)+K-1])
    alpha <- with(env, threshold(theta,K))
    env$alpha <- alpha
    env$beta <- beta
    if (env$p>0) eta <- env$X%*%beta else eta <- cbind(rep(0,env$n))
    env$lp <- kronecker(-eta,rbind(alpha),"+")
    F <- with(env,h(lp))
    if (tolower(type)[1]=="cumulative") return(F)            
    Pr <- cbind(F,1)-cbind(0,F)
    colnames(Pr) <- object$levels
    return(Pr)
}
