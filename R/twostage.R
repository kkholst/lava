twostagelvm <- function(object, model2,
                formula=NULL, model.object=FALSE, predict.fun=NULL,
                type="quadratic",...) {
    if (!inherits(model2,c("lvm"))) stop("Expected lava object ('lvm',...)")
    if (!is.null(formula)) {
        model2 <- nonlinear(model2, formula, type=type)
    }
    nonlin <- NULL
    val <- nonlinear(model2)
    if (is.null(formula) && length(val)==0 && length(nonlinear(object))>0) {
        val <- nonlinear(object)
    }
    xnam <- c()
    if (length(val)>0) {
        predict.fun <- NULL
        for (i in seq_along(val)) {
            if (!all(val[[i]]$newx%in%xnam)) {
                xnam <- union(xnam,val[[i]]$newx)
                predict.fun <- c(predict.fun, list(val[[i]]$pred))
            }
            model2$attributes$nonlinear <- NULL
            if (inherits(object,"lvmfit")) {
                object$model$attributes$nonlinear <- NULL
            }
            model2 <- regression(model2, to=names(val)[i], from=val[[i]]$newx)
        }
        nonlin <- val
    }
    if (model.object) {
        model <- Model(object) %++% model2
        cl <- match.call(expand.dots=TRUE)
        cl[[1]] <- twostage
        cl$object <- object
        cl$model2 <- model2
        cl$predict.fun <- predict.fun
        cl["model.object"] <- NULL
        return(structure(list(model=model, nonlinear=nonlin, call=cl), class="twostage.lvm"))
    }
    res <- c(list(object=object, model2=model2), list(...))
    res$predict.fun <- predict.fun
    res$nonlinear <- val
    return(res)
}


uhat <- function(p=coef(model1), model1, data=model.frame(model1), nlobj) {
    if (!is.function(nlobj)) {
        predict.fun <- lapply(nlobj, function(x) x[["pred"]])
    } else { predict.fun <- nlobj }
    if (inherits(model1, "lvm.mixture")) {
        if (is.list(predict.fun)) {
            unams <- lapply(nlobj,function(x) x$newx)
            unam <- unique(unlist(unams))
            res <- matrix(0, NROW(data), ncol=length(unam))
            colnames(res) <- unam
            for (i in seq_along(predict.fun)) {
                res[, unams[[i]]] <-
                    predict(model1, p=p, data=data, predict.fun=predict.fun[[i]])
            }
            return(res)
        } else {
            Pr <- cbind(predict(model1, p=p, data=data, predict.fun=predict.fun))
            return(Pr)
        }
        ##P <- list(mean=Pr, var=attr(Pr,"cond.var"))
    }  else {
        P <- predictlvm(model1, p=p, data=data)
    }
    if (is.list(predict.fun)) {
        unams <- lapply(nlobj,function(x) x$newx)
        unam <- unique(unlist(unams))
        args <- list(P$mean, P$var, data)
        res <- matrix(0, NROW(data), ncol=length(unam))
        colnames(res) <- unam
        for (i in seq_along(predict.fun)) {
            res[, unams[[i]]] <- do.call(predict.fun[[i]], args)
        }
        return(res)
    }
    return(cbind(predict.fun(P$mean, P$var, model.frame(model1))))
}


##' Two-stage estimator
##'
##' Generic function.
##'
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
##' @param predict.fun Prediction of latent variable
##' @param id1 Optional id-variable (stage 1 model)
##' @param id2 Optional id-variable (stage 2 model)
##' @param all If TRUE return additional output (naive estimates)
##' @param formula optional formula specifying non-linear relation
##' @param std.err If FALSE calculations of standard errors will be skipped
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
##' d$u2 <- pp1[,]^2+attr(pp1,"cond.var")[1]
##' m2 <- lvm(c(y1,y2,y3)~eta,c(y1,eta)~u1+u2+z); latent(m2) <- ~eta
##' e2 <- estimate(m2,d)
##' }
##'
##' ## Two-stage
##' m1 <- lvm(c(x1,x2,x3)~f1,f1~z); latent(m1) <- ~f1
##' m2 <- lvm(c(y1,y2,y3)~eta,c(y1,eta)~u1+u2+z); latent(m2) <- ~eta
##' pred <- function(mu,var,data,...)
##'     cbind("u1"=mu[,1],"u2"=mu[,1]^2+var[1])
##' (mm <- twostage(m1,model2=m2,data=d,predict.fun=pred))
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
##' nonlinear(m2,type="spline",knots=kn) <- eta2~eta1
##' a <- twostage(m1,m2,data=d)
##' x <- seq(-3,3,by=0.1)
##' y <- predict(a, newdata=data.frame(eta1=x))
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
##'
##' \dontrun{ ## Reduce timing
##'  ## Cross-validation example
##'  ma <- lvm(c(x1,x2,x3)~u,latent=~u)
##'  ms <- functional(ma, y~u, f=function(x) -.4*x^2)
##'  d <- sim(ms,500)#,seed=1)
##'  ea <- estimate(ma,d)
##'
##'  mb <- lvm()
##'  mb1 <- nonlinear(mb,type="linear",y~u)
##'  mb2 <- nonlinear(mb,type="quadratic",y~u)
##'  mb3 <- nonlinear(mb,type="spline",knots=c(-3,-1,0,1,3),y~u)
##'  mb4 <- nonlinear(mb,type="spline",knots=c(-3,-2,-1,0,1,2,3),y~u)
##'  ff <- lapply(list(mb1,mb2,mb3,mb4),
##'	      function(m) function(data,...) twostage(ma,m,data=data,st.derr=FALSE))
##'  a <- cv(ff,data=d,rep=1,mc.cores=1)
##'  a
##'}
twostage.lvmfit <- function(object, model2, data=NULL,
                    predict.fun=NULL,
                    id1=NULL, id2=NULL, all=FALSE,
                    formula=NULL, std.err=TRUE,
                    ...) {
    if (!is.null(predict.fun)) {
        object$attributes$nonlinear <- list()
        model2$attributes$nonlinear <- list()
    }
    val <- twostagelvm(object=object,model2=model2,predict.fun=predict.fun,
                      id1=id1, id2=id2, all=all, formula=formula, ...)
    object <- val$object
    model2 <- val$model2
    predict.fun <- val$predict.fun
    p1 <- coef(object)
    if (length(val$nonlinear)==0) {
        val$nonlinear <- predict.fun
    }
    pp <- uhat(p1,object,nlobj=val$nonlinear)
    newd <- data
    newd[,colnames(pp)] <- pp

    model2 <- estimate(model2,data=newd,...)
    p2 <- coef(model2)
    if (std.err) {
        if (is.null(id1)) id1 <- seq(nrow(model.frame(object)))
        if (is.null(id2)) id2 <- seq(nrow(model.frame(model2)))
        model1 <- object
        if (!inherits(object,"estimate")) {
            model1 <- estimate(NULL,coef=p1,id=id1,iid=iid(object))
        }

        e2 <- estimate(model2, id=id2)
        U <- function(alpha=p1,beta=p2) {
            pp <- uhat(alpha,object,nlobj=val$nonlinear)
            newd <- model.frame(model2)
            newd[,colnames(pp)] <- pp
            score(model2,p=beta,data=newd)
        }
        Ia <- numDeriv::jacobian(function(p) U(p),p1)
        stacked <- stack(model1,e2,Ia)
    } else {
        e2 <- estimate(coef=p2,vcov=NA)
    }
    coef <- model2$coef
    res <- model2
    res$estimator <- "generic"

    if (std.err) {
        res[names(stacked)] <- stacked
        cc <- stacked$coefmat[,c(1,2)];
        cc <- cbind(cc,cc[,1]/cc[,2],stacked$coefmat[,5])
        coef[,] <- cc
        res$coef <- coef
        res$vcov <- vcov(stacked)
        if (all) {
            res$naive <- model2
            res$naive.robust <- e2
        }
    } else {
        res$coef[,-1] <- NA
    }
    res$fun <- predict.fun
    res$estimate1 <- object
    res$estimate2 <- model2
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
    m2$attributes$nonlinear <- nl
    twostage(object=m1,model2=m2,data=data,predict.fun=nl[[1]]$pred,...)
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
predict.twostage.lvmfit <- function(object,
                            newdata,
                            variable=names(nonlinear(object)),
                            p=coef(object),
                            type=c("model2","latent"),
                            x=NULL,
                            ...) {
    if (missing(newdata)) stop("provide data for prediction")
    nl <- nonlinear(object)
    unam <- unique(unlist(lapply(nl,function(x) x$x)))
    if (is.vector(newdata) || all(unam%in%colnames(newdata)))
        type <- "latent"


    if (tolower(type[1])%ni%c("latent")) {
        p1 <- coef(object$estimate1)
        pred1 <- uhat(p1, data=newdata, object$estimate1, nlobj=nl)
        if (tolower(type[1])==c("model1"))
            return(pred1)
        newdata <- as.data.frame(newdata)
        newdata[,colnames(pred1)] <- pred1
        pred <- predict(object$estimate2,x=x,p=p,data=newdata,...)
        attr(pred,"p") <- NULL
        attr(pred,"e") <- NULL
        attr(pred,"cond.var") <- NULL
        return(pred)
    }

    ## Association between predicted latent variables and child nodes:
    if (is.numeric(variable)) {
        variable <- names(nonlinear(object))[variable]
    }
    nl <- nl[variable]
    res <- matrix(nrow=NROW(newdata),ncol=length(nl))
    colnames(res) <- names(nl)
    ##unam <- unique(unlist(lapply(nl, function(x) x$newx)))
    ##colnames(res) <- unam
    if (!is.null(x)) {
        newd <- newdata
        m0 <- Model(object$estimate2)
        p0 <- p
    }
    for (i in seq_along(nl)) {
        pnam <- c(variable,paste0(variable,"~",nl[[i]]$newx))
        pidx <- match(pnam,names(coef(object)))
        b <- p[pidx]
        F <- nl[[i]]$f
        if (is.vector(newdata)) {
            res[,i] <- F(b,newdata)
        } else {
            res[,i] <- F(b,newdata[,nl[[i]]$x])
            if (!is.null(x)) {
                newd[,nl[[i]]$newx] <- 0
                latent(m0,clear=TRUE) <- names(nl)[i]
                regression(m0, y=names(nl)[i], x=paste0(names(nl)[i],".offset")) <- 1
                p0[pidx] <- 0
            }
        }
    }

    if (!is.null(x) && !is.vector(newdata)) {
        latentres <- res
        colnames(latentres) <- paste0(colnames(res),".offset",collapse="")
        newd <- cbind(newd,latentres)
        res <- predict(m0,data=newd,p=p0,x=x,...)
        attr(res,"cond.var") <- NULL

    }
    return(res)
}

##' Cross-validated two-stage estimator
##'
##' Cross-validated two-stage estimator for non-linear SEM
##' @export
##' @param model1 model 1 (exposure measurement error model)
##' @param model2 model 2
##' @param data data.frame
##' @param control1 optimization parameters for model 1
##' @param control2 optimization parameters for model 1
##' @param knots.boundary boundary points for natural cubic spline basis
##' @param mc.cores number of cores to use for parallel computations
##' @param nmix number of mixture components
##' @param df spline degrees of freedom
##' @param fix automatically fix parameters for identification (TRUE)
##' @param std.err calculation of standard errors (TRUE)
##' @param nfolds Number of folds (cross-validation)
##' @param rep Number of repeats of cross-validation
##' @param messages print information (>0)
##' @param ... additional arguments to lower level functions
##' @examples
##' \donttest{ ## Reduce Ex.Timings
##' m1 <- lvm( x1+x2+x3 ~ u1, latent= ~u1)
##' m2 <- lvm( y ~ 1 )
##' m <- functional(merge(m1,m2), y ~ u, f=function(x) sin(x)+x)
##' distribution(m, ~u1) <- uniform.lvm(-6,6)
##' d <- sim(m,n=500,seed=1)
##' nonlinear(m2) <- y~u1
##' val <- twostageCV(m1, m2, data=d, std.err=FALSE, df=2:6, nmix=1:2,
##'                   nfolds=2, mc.cores=1)
##' val
##' }
twostageCV <- function(model1, model2, data, control1=list(trace=0), control2=list(trace=0),
                knots.boundary, mc.cores=1, nmix=1:4, df=1:9, fix=TRUE, std.err=TRUE,
                nfolds=5, rep=1, messages=0, ...) {
    op <- options(warn=-1)
    if (fix) {
        model1 <- baptize(fixsome(model1, param="relative"))
        intercept(model1, latent(model1)) <- NA
    }
    e1a <- estimate(model1, data=data, control=control1)
    if (missing(knots.boundary))
        knots.boundary <- range(predict(e1a,vars(e1a)))
    ## Starting values for mixture models
    plab <- parlabels(model1)
    pfree <- setdiff(coef(model1),plab)
    pfree.idx <- match(pfree,coef(model1)) ## Index of unlabeled parameters
    intpos <- setdiff(parpos(model1)$v,0)
    pfree.int <- intersect(pfree.idx,intpos) ## Free intercept parameters
    p0 <- coef(e1a)
    pint <- p0[setdiff(intpos,pfree.int)]
    startf <- function(n) {
        u0 <- seq(knots.boundary[1],knots.boundary[2],length.out=n); names(u0) <- paste0("p",seq_along(u0))
        c(pint,u0,p0[-intpos])
    }
    control1$start <- NULL
    ee <- list(e1a)
    nmix <- setdiff(nmix,1)
    for (k in nmix) {
        if (messages>0) cat("Fitting mixture model with", k, "components\n")
        ee <- c(ee, list(mixture(model1, k=k, data=data,
                           control=c(control1,list(start=startf(k))))))
    }
    AIC1 <- unlist(lapply(ee,AIC))
    names(AIC1) <- c(1,nmix)
    ii <- which.min(AIC1)
    ## Exposure measurement model
    e1 <- ee[[ii]] ## Selected model by AIC

    F <- nonlinear(model2)
    form <- as.formula(paste(names(F)[1], "~", x=F[[1]]$x))

    MM <- list(nonlinear(model2, form, type="linear"))
    df <- setdiff(df, 1)
    Knots <- list()
    for (i in df) {
        knots <- seq(knots.boundary[1],knots.boundary[2],length.out=i+1)
        Knots  <- c(Knots, list(knots))
        MM <- c(MM,
               list(nonlinear(model2, form, type="spline", knots=knots)))
    }
    if (!inherits(e1, "lvm.mixture")) {
        f0 <- function(data) list(e0=estimate(model1,data=data,control=c(control1,list(start=coef(e1)))))
    } else {
        f0 <- function(data) list(e0=mixture(model1,data=data,k=e1$k,control=c(control1,list(start=coef(e1)))))
    }
    ff <- lapply(MM,
                 function(m) function(data,e0,...)  twostage(e0,m,data=data,std.derr=FALSE))
    a <- cv(ff,data=data,K=nfolds,rep=rep,mc.cores=mc.cores,shared=f0)
    M <- MM[[which.min(coef(a))]]
    e2 <- twostage(e1,M,data,control=control2, std.err=std.err)
    options(op)
    res <- list(AIC1=cbind(AIC1), model1=e1, model2=e2, cv=coef(a), knots=c(list(NA),Knots),
                nfolds=nfolds, rep=rep)
    structure(res, class="twostageCV")
}


##' @export
print.twostageCV <- function(x,...) {
    printline(70)
    i1 <- which.min(x$AIC1)
    nmix <- rownames(x$AIC1)[i1]
    cat("Selected mixture model: ",nmix," component", ifelse(i1>1, "s",""),"\n", sep="")
    print(x$AIC1)
    printline(70)
    i2 <- which.min(x$cv)
    splinedf <- unlist(lapply(x$knots,function(x) if (any(is.na(x))) return(1) else length(x)-1))
    cat("Selected spline model degrees of freedom: ", splinedf[i2] ,"\n", sep="")
    knots <- rbind(x$knots[[i2]])
    cat("Knots:", paste(formatC(knots,...) , collapse=" "), "\n\n")
    rmse <- x$cv
    rownames(rmse) <- paste0("df:",splinedf)
    colnames(rmse) <- paste0("RMSE(nfolds=",x$nfolds,", rep=",x$rep,")")
    print(rmse)
    printline(70)
    cat("\n")
    print(CoefMat(x$model2,...),quote=FALSE)
}

##' @export
coef.twostageCV <- function(object,...) {
    coef(Model(object),...)
}

vcov.twostageCV <- function(object,...) {
    vcov(Model(object),...)
}

iid.twostageCV <- function(object,...) {
    iid(Model(object),...)
}

Model.twostageCV <- function(x,...) {
    x$model2
}


##' @export
summary.twostageCV <- function(object,...) {
    with(object, list(model1=summary(model1),
                      AIC1=AIC1, cv=cv, knots=knots,
                      model2=summary(model2)))
}

##' @export
predict.twostageCV <- function(object,... ) {
    predict(Model(object),...)
}
