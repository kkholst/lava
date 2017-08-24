twostagelvm <- function(object, model2, 
                formula=NULL, model.object=FALSE, predictfun=NULL,
                type="quadratic",...) {
    if (!inherits(model2,c("lvm"))) stop("Expected lava object ('lvm',...)")
    if (!is.null(formula)) {
        model2 <- nonlinear(model2, formula, type=type)
    }
    nonlin <- NULL
    val <- nonlinear(model2)
    if (length(nonlinear(object))>0) {
        val <- nonlinear(object)
    }
    if (length(val)>0) {
        predictfun <- val[[1]]$pred
        model2$attributes$nonlinear <- NULL
        if (inherits(object,"lvmfit")) {
            object$model$attributes$nonlinear <- NULL
        }
        model2 <- regression(model2, to=names(val)[1], from=val[[1]]$newx)
        nonlin <- val[1]
    }
    if (model.object) {
        model <- Model(object) %++% model2
        cl <- match.call(expand.dots=TRUE)
        cl[[1]] <- twostage
        cl$object <- object
        cl$model2 <- model2
        cl$predictfun <- predictfun
        cl["model.object"] <- NULL
        return(structure(list(model=model, nonlinear=nonlin, call=cl), class="twostage.lvm"))
    }
    res <- c(list(object=object, model2=model2), list(...))
    res$predictfun <- predictfun
    res$nonlinear <- val
    return(res)
}


##' Two-stage estimator
##'
##' Generic function
##' @seealso twostage.lvm twostage.lvmfit twostage.lvm.mixture twostage.estimate
##' @export
##' @param object Model object
##' @param ... Additional arguments to lower level functions
"twostage" <- function(object,...) UseMethod("twostage")

##' Two-stage estimator (non-linear SEM)
##'
##' Two-stage estimator for non-linear structural equation models
##' @export
##' @param object Stage 1 measurement model
##' @param model2 Stage 2 SEM
##' @param data data.frame
##' @param predictfun Prediction of latent variable
##' @param id1 Optional id-variable (stage 1 model)
##' @param id2 Optional id-variable (stage 2 model)
##' @param all If TRUE return additional output (naive estimates)
##' @param formula optional formula specifying non-linear relation
##' @param ... Additional arguments to lower level functions
##' @aliases twostage.lvmfit twostage.lvm twostage.lvm.mixture twostage.estimate nonlinear nonlinear<-
##' @examples
##' m <- lvm(c(x1,x2,x3)~f1,f1~z,
##'          c(y1,y2,y3)~f2,f2~f1+z)
##' latent(m) <- ~f1+f2
##' d <- simulate(m,100,p=c("f2,f2"=2,"f1,f1"=0.5),seed=1)
##'
##' ## Full MLE
##' ee <- estimate(m,d)
##'
##' ## Manual two-stage
##' \dontrun{
##' m1 <- lvm(c(x1,x2,x3)~f1,f1~z); latent(m1) <- ~f1
##' e1 <- estimate(m1,d)
##' pp1 <- predict(e1,f1~x1+x2+x3)
##'
##' d$u1 <- pp1[,]
##' d$u2 <- pp1[,]^2+attr(pp1,"cond.var")
##' m2 <- lvm(c(y1,y2,y3)~eta,c(y1,eta)~u1+u2+z); latent(m2) <- ~eta
##' e2 <- estimate(m2,d)
##' }
##'
##' ## Two-stage
##' m1 <- lvm(c(x1,x2,x3)~f1,f1~z); latent(m1) <- ~f1
##' m2 <- lvm(c(y1,y2,y3)~eta,c(y1,eta)~u1+u2+z); latent(m2) <- ~eta
##' pred <- function(mu,var,data,...)
##'     cbind("u1"=mu[,1],"u2"=mu[,1]^2+var[1])
##' (mm <- twostage(m1,model2=m2,data=d,predictfun=pred))
##'
##' if (interactive()) {
##'     pf <- function(p) p["eta"]+p["eta~u1"]*u + p["eta~u2"]*u^2
##'     plot(mm,f=pf,data=data.frame(u=seq(-2,2,length.out=100)),lwd=2)
##' }
##'
##' ## Splines
##' f <- function(x) cos(2*x)+x+-0.25*x^2
##' m <- lvm(x1+x2+x3~eta1, y1+y2+y3~eta2, latent=~eta1+eta2)
##' functional(m, eta2~eta1) <- f
##' d <- sim(m,500,seed=1,latent=TRUE)
##' m1 <- lvm(x1+x2+x3~eta1,latent=~eta1)
##' m2 <- lvm(y1+y2+y3~eta2,latent=~eta2)
##' mm <- twostage(m1,m2,formula=eta2~eta1,type="spline")
##' if (interactive()) plot(mm)
##'
##' nonlinear(m2,type="quadratic") <- eta2~eta1
##' a <- twostage(m1,m2,data=d)
##' if (interactive()) plot(a)
##'
##' kn <- c(-1,0,1)
##' nonlinear(m2,type="spline",knots=kn,boundary=c(-4,4)) <- eta2~eta1
##' a <- twostage(m1,m2,data=d)
##' x <- seq(-3,3,by=0.1)
##' y <- predict(a, newdata=x)
##'
##' if (interactive()) {
##'   plot(eta2~eta1, data=d)
##'   lines(x,y, col="red", lwd=5)
##'
##'   p <- estimate(a,f=function(p) predict(a,p=p,newdata=x))$coefmat
##'   plot(eta2~eta1, data=d)
##'   lines(x,p[,1], col="red", lwd=5)
##'   confband(x,lower=p[,3],upper=p[,4],center=p[,1], polygon=TRUE, col=Col(2,0.2))
##'
##'   l1 <- lm(eta2~splines::ns(eta1,knots=kn),data=d)
##'   p1 <- predict(l1,newdata=data.frame(eta1=x),interval="confidence")
##'   lines(x,p1[,1],col="green",lwd=5)
##'   confband(x,lower=p1[,2],upper=p1[,3],center=p1[,1], polygon=TRUE, col=Col(3,0.2))
##' }
twostage.lvmfit <- function(object, model2, data=NULL,
                    predictfun=function(mu,var,data,...)
                        cbind("u1"=mu[,1],"u2"=mu[,1]^2+var[1]),
                    id1=NULL, id2=NULL, all=FALSE,
                    formula=NULL, ...) {
    val <- twostagelvm(object=object,model2=model2,predictfun=predictfun,
                      id1=id1, id2=id2, all=all, formula=formula, ...)
    object <- val$object
    model2 <- val$model2
    predictfun <- val$predictfun
    p1 <- coef(object)
    uhat <- function(p=p1) {
        if (inherits(object,"lvm.mixture")) {
            Pr <- predict(object,p=p)
            P <- list(mean=Pr,var=attr(Pr,"cond.var"))
        }  else {
            P <- predictlvm(object,p=p,data=model.frame(object))
        }
        return(cbind(predictfun(P$mean,P$var,model.frame(object))))
    }

    pp <- uhat()
    newd <- data
    newd[,colnames(pp)] <- pp
    model2 <- estimate(model2,data=newd,...)
    p2 <- coef(model2)
    if (is.null(id1)) id1 <- seq(nrow(model.frame(object)))
    if (is.null(id2)) id2 <- seq(nrow(model.frame(model2)))
    model1 <- object
    if (!inherits(object,"estimate")) {
        model1 <- estimate(NULL,coef=p1,id=id1,iid=iid(object))
    }
    e2 <- estimate(model2, id=id2)
    U <- function(alpha=p1,beta=p2) {
        pp <- uhat(alpha)
        newd <- model.frame(model2)
        newd[,colnames(pp)] <- pp
        score(model2,p=beta,data=newd)
    }
    Ia <- -numDeriv::jacobian(function(p) U(p),p1)
    stacked <- stack(model1,e2,Ia)

    coef <- model2$coef
    res <- model2
    res[names(stacked)] <- stacked
    res$estimator <- "generic"
    cc <- stacked$coefmat[,c(1,2)];
    cc <- cbind(cc,cc[,1]/cc[,2],stacked$coefmat[,5])
    coef[,] <- cc
    res$coef <- coef
    res$fun <- predictfun
    if (all) {
        res$naive <- model2
        res$naive.robust <- e2
    }
    res$nonlinear <- val$nonlinear
    structure(res,class=c("twostage.lvmfit","measurement.error","lvmfit","estimate"))
}

##' @export
estimate.twostage.lvm <- function(x,data,...) {
    if (missing(data)) stop("'data' needed")
    m1 <- x$call$object
    m2 <- x$call$model2
    nl <- x$nonlinear
    if (!inherits(m1,"lvmfit")) {
        args <- c(list(x=m1, data=data), list(...))
        args <- args[intersect(names(as.list(base::args(estimate.lvm))),names(args))]
        m1 <- do.call(estimate, args)
    }
    twostage(object=m1,model2=m2,data=data,predictfun=nl[[1]]$pred,...)
}

##' @export
twostage.twostage.lvm <- function(object,...) estimate.twostage.lvm(object,...)


##' @export
twostage.lvm <- function(object,model2,data=NULL, ...) {
    if (is.null(data)) {
        return(twostagelvm(object=object, model2=model2, model.object=TRUE, ...))
    }
    args <- c(list(x=object, data=data), list(...))
    args <- args[intersect(names(as.list(base::args(estimate.lvm))),names(args))]
    e1 <- do.call(estimate, args)
    twostage(object=e1,model2=model2,data=data, ...)    
}

##' @export 
twostage.lvm.mixture <- twostage.lvmfit

##' @export
twostage.estimate <- twostage.lvmfit

##' @export
print.twostage.lvm <- function(x,...) {
    printline()
    cat("Model 1:\n")
    print(Model(x$call$object))
    printline()
    cat("Model 2:\n")
    print(Model(x$call$model2))

}

##' @export
plot.twostage.lvm <- function(x,...) {
    model <- x$model
    m1 <- Model(x$call$object)
    m2 <- x$call$model2
    nl <- nonlinear(x)
    model <- regression(model, to=nl[[1]]$newx, from=nl[[1]]$x)
    elist <- edgeList(m1)
    vlist <- vars(m1)
    model <- beautify(model)
    for (i in seq_len(nrow(elist))) {
        e <- toformula(y=vlist[elist[i,2]],x=vlist[elist[i,1]])
        edgelabels(model, e, cex=0.7) <- 1
    }
    elist <- edgeList(m2)
    vlist <- vars(m2)
    for (i in seq_len(nrow(elist))) {
        e <- toformula(vlist[elist[i,2]],vlist[elist[i,1]])
        edgelabels(model, e, cex=0.7) <- 2
    }
    nodecolor(model, nl[[1]]$newx) <- "gray"
    for (xx in nl[[1]]$newx) {
        e <- toformula(y=names(nl)[1],x=xx)
        edgelabels(model,e,col="gray", cex=0.7, lty=1) <- 2
    }
    for (xx in nl[[1]]$newx) {
        e <- toformula(y=xx,x=nl[[1]]$x)
        edgelabels(model,e,col="gray", cex=0.7, lty=2) <- ""
    }
    plot(model, ...)    
}


##' @export
predict.twostage.lvmfit <- function(object,p=coef(object),variable=1,newdata,...) {
    if (missing(newdata)) stop("provide data for prediction")
    if (is.numeric(variable)) {
        variable <- names(nonlinear(object))[variable]
    }
    nl <- nonlinear(object)[[variable]]
    pnam <- c(variable,paste0(variable,"~",nl$newx))
    pidx <- match(pnam,names(coef(object)))
    b <- p[pidx]
    F <- nl$f
    yhat <- F(b,newdata)
    return(yhat)
}
