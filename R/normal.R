normal_method.lvm <- "nlminb0"

normal_objective.lvm <- function(x,p,data,weight2=NULL,indiv=FALSE,...) {
    save.seed <- .Random.seed; set.seed(1)
    on.exit(set.seed(save.seed))
    y.idx <- lava::index(x)$endo.idx
    mu <- predict(x,data=data,p=p)
    S <- attributes(mu)$cond.var
    class(mu) <- "matrix"
    y <- lava::endogenous(x)
    ord <- lava::ordinal(x)
    status <- rep(0,length(y))
    if (exists("binary.lvm")) status[match(do.call("binary",x),y)] <- 2
    status[match(ord,y)] <- 2
    thres <- matrix(0,nrow=length(y),max(1,attributes(ord)$K-1)); rownames(thres) <- y
##    browser()
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
    if (!is.null(weight2)) {
        yu[,colnames(weight2)] <- weight2
        status[match(colnames(weight2),y)] <- 1
    }

    l <- mets::loglikMVN(yl,yu,status,mu,S,thres)
    ## l <- .Call("loglikMVN",
    ##            yl=as.matrix(yl),
    ##            yu=as.matrix(yu),
    ##            status=as.integer(status),
    ##            mu=as.matrix(mu),dmu=NULL,s=as.matrix(S),ds=NULL,
    ##            z=NULL,su=NULL,dsu=NULL,
    ##            threshold=as.matrix(thres),
    ##            dthreshold=NULL, package="mets")
  
    if (indiv) return(-l)
    return(-sum(l))  
}

normal_logLik.lvm <- function(object,p,data,weight2,...) {
    res <- -normal_objective.lvm(x=object,p=p,data=data,weight2=weight2,...)
    return(res)
}

normal_gradient.lvm <- function(x,p,data,weight2,indiv=FALSE,...) {
    if (indiv) {
        return(numDeriv::jacobian(function(p0) normal_objective.lvm(x,p=p0,data=data,weight2=weight2,indiv=TRUE,...),p,method=lava.options()$Dmethod))
    }
    numDeriv::grad(function(p0) normal_objective.lvm(x,p=p0,data=data,weight2=weight2,...),p,method=lava.options()$Dmethod)
}

