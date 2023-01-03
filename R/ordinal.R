ordinal.remove.hook <- function(x,var,...) {
    ordinal(x,K=0) <- var
    return(x)
}

color.ordinal <- function(x,subset=vars(x),...) {
    return(list(vars=intersect(subset,ordinal(x)),col="indianred1"))
}

ordinal.sim.hook <- function(x,data,p,modelpar,...) {
    ovar <- ordinal(x)

    for (i in seq_len(length(ovar))) {
        if (attributes(ovar)$liability[i]) {
            idx <- attributes(ovar)$idx[[ovar[i]]]
            if (length(idx)==0) {
                breaks <- c(-Inf,0,Inf)
            } else {
                breaks <- c(-Inf,ordreg_threshold(modelpar$e[idx]),Inf)
            }
            z <- cut(data[,ovar[i]],breaks=breaks)
            data[,ovar[i]] <- as.numeric(z)-1
        }
        K <- attributes(ovar)$K[i]
        lab <- attributes(ovar)$labels[ovar[i]][[1]]
        if (!is.null(lab))
            data[,ovar[i]] <- factor(data[,ovar[i]],
                                     levels=seq(K)-1,
                                     labels=lab)
    }
    return(data)
}

ordinal.estimate.hook <- function(x,data,weights,data2,estimator,...) {
    dots <- list(...)

    nestimator <- c("normal")
    nestimator2 <- c("tobit","tobitw","gaussian")

    ord <- ordinal(x)
    bin <- NULL
    if (is.null(estimator) && length(ord)>0) estimator <- nestimator[1]

    ## Binary outcomes -> censored regression
    if (is.null(dim(data))) return(NULL)
    if (is.null(estimator) || estimator%in%c(nestimator2,nestimator)) {
        for (i in setdiff(lava::endogenous(x),bin)) {
            if (is.character(data[,i]) | is.factor(data[,i])) { # Transform binary 'factor'
                y <- as.factor(data[,i])
                data[,i] <- as.numeric(y)-1
                estimator <- nestimator[1]
                ordinal(x,K=nlevels(y)) <- i
            }
        }
        ord <- ordinal(x)
        if (length(ord)>0 && !is.null(estimator) && estimator%in%nestimator2) {
          estimator <- nestimator[1]
        }
    }

    ## Transform 'Surv' objects
    data2 <- mynames <- NULL
    if (is.null(estimator) || estimator%in%nestimator[1] ) {
        for (i in setdiff(lava::endogenous(x),c(bin,ord))) {
            if (survival::is.Surv(data[,i])) {
                S <- data[,i]
                y1 <- S[,1]
                if (attributes(S)$type=="left")  {
                    y2 <- y1
                    y1[S[,2]==0] <- -Inf
                }
                if (attributes(S)$type=="right") {
                    y2 <- y1
                    y2[S[,2]==0] <- Inf
                }
                if (attributes(S)$type=="interval2") {
                    y2 <- S[,2]
                }
                if (attributes(S)$type=="interval") {
                    y2 <- S[,2]
                    y2[S[,3]==1L] <- y1[S[,3]==1L]
                }
                if (!(attributes(S)$type%in%c("left","right","interval2","interval"))) stop("Surv type not supported.")
                mynames <- c(mynames,i)
                y2 <- cbind(y2)
                colnames(y2) <- i
                data2 <- cbind(data2,y2)
                data[,i] <- y1
                estimator <- "normal"
            }
        }
    }

    return(c(list(x=x,data=data,weights=weights,data2=data2,estimator=estimator),dots))
}


##' Define variables as ordinal
##'
##' Define variables as ordinal in latent variable model object
##' @export
##' @aliases ordinal ordinal<-
##' @param x Object
##' @param ... additional arguments to lower level functions
##' @param value variable (formula or character vector)
##' @examples
##' if (requireNamespace("mets")) {
##' m <- lvm(y + z ~ x + 1*u[0], latent=~u)
##' ordinal(m, K=3) <- ~y+z
##' d <- sim(m, 100, seed=1)
##' e <- estimate(m, d)
##' }
##' 
"ordinal<-" <- function(x,...,value) UseMethod("ordinal<-")

##' @export
"ordinal<-.lvm" <- function(x,...,value) {
  ordinal(x, value, ...)
}

##' @export
"ordinal" <- function(x,...) UseMethod("ordinal")

##' @export
print.ordinal.lvm <- function(x,...) {
  cat(rep("_",28),"\n",sep="")
  for (i in x) {
    val <- attr(x,"fix")[[i]]
    if (length(val)==0)
      cat(paste(i,"binary",sep=":"),"\n")
    else print(unlist(attr(x,"fix")[[i]]),quote=FALSE)
    cat(rep("_",28),"\n",sep="")
  }
}

##' @export
`ordinal.lvm` <- function(x,var=NULL,K=2, constrain, breaks=NULL, p, liability=TRUE, labels, exo=FALSE, ...) {
    if (inherits(var,"formula")) {
        var <- all.vars(var)
    }
    if (is.null(var)) {
        ordidx <- unlist(x$attributes$ordinal)
        KK <- unlist(x$attributes$nordinal)
        idx <- x$attributes$ordinalparname
        fix <- lapply(idx,function(z) x$exfix[z])
        liability <- x$attributes$liability
        labels <- x$attributes$labels
        if (length(ordidx)>0) {
            val <- names(ordidx)
            return(structure(val,K=KK,idx=idx,fix=fix,liability=liability,labels=labels,class="ordinal.lvm"))
        }
        else
            return(NULL)
    }
    if (K[1]==0L || is.null(K[1]) || (is.logical(K) & !K[1])) {
        idx <- na.omit(match(var,names(x$attributes$ordinal)))
        if (length(idx)>0) {
            pp <- unlist(x$attributes$ordinalparname[idx])
            if (!is.null(pp)) parameter(x,remove=TRUE) <- pp
            if (!is.null(x$attributes$ordinalparname))
                x$attributes$ordinalparname <- x$attributes$ordinalparname[-idx]
            x$attributes$ordinal <- x$attributes$ordinal[-idx]
            ##x$attributes$labels[var] <- NULL
            x$attributes$type <- x$attributes$type[-idx]
            x$attributes$constrainY <- x$attributes$constrainY[setdiff(names(x$attributes$constrainY),var)]
            x$attributes$liability <- x$attributes$liability[-idx]
            x$attributes$nordinal <- x$attributes$nordinal[-idx]
            x$attributes$normal <- x$attributes$normal[-idx]
            exo <- intersect(var,exogenous(x,latent=TRUE))
            if (length(exo)>0) {
                intercept(x,var) <- NA
                covariance(x,var) <- NA
                exogenous(x) <- union(exogenous(x),exo)
            }
        }
        return(x)
    }
    
    if (!missing(p)) breaks <- qnorm(cumsum(p))
    if (!is.null(breaks)) {
        breaks <- ordreg_ithreshold(breaks)
        K <- length(breaks)+1
    }
    if (!missing(labels)) K <- length(labels)
    if (length(var)>length(K)) K <- rep(K[1],length(var))
    if (length(var)==1 && !missing(constrain)) constrain <- list(constrain)
    if (length(var)>1) {
        if (!missing(labels) && !is.list(labels)) labels <- rep(list(labels),length(var))
        if (!missing(breaks) && !is.list(breaks)) breaks <- rep(list(breaks),length(var))
        if (!missing(constrain) && !is.list(constrain)) constrain <- rep(list(constrain),length(var))
    }

    addvar(x) <- var
    for (i in seq_len(length(var))) {
        if (K[i]>2 || (K[i]==2 && !liability)) {
            parname <- paste0(var[i],":",paste(seq(K[i]-1)-1,seq(K[i]-1),sep="|"))
            newpar <- if (is.null(breaks)) {
                rep(-1,K[i]-1)
            } else if (is.list(breaks)) breaks[[i]] else breaks
            if (length(newpar)<K[i]-1) stop("Wrong number of starting values")
            newfix <- if (missing(constrain))
                rep(list(NA),length(newpar)) else constrain[[i]]
            if (length(newfix)<K[i]-1) stop("Wrong number of constraints")
            names(newpar) <- names(newfix) <- parname

            parameter(x,newfix,start=newpar) <- names(newfix)
            ## pex <- parname%in%names(x$expar)
            ## if (any(pex <- parname%in%names(x$expar))) {
            ##     if (!all(pex)) stop("Cannot change number of categories! Re-specify model.")
            ##     x$attributes$iordinal[var] <- list(idx)
            ## }
            x$attributes$ordinalparname[var[i]] <- list(names(newfix))
        }
        x$attributes$type[var[i]] <- ifelse(K[i]>2,"categorical","binary")
        if (K[i]>2) intfix(x,var[i],NULL) <- 0
        if (!liability) {
            mytr <- function(y,p,idx,...) {
                breaks <- c(-Inf,ordreg_threshold(p[idx]),Inf)
                as.numeric(cut(y,breaks=breaks))-1
            }
            myalist <- substitute(alist(y=,p=,idx=pp),
                                  list(pp=x$attributes$ordinalparname[[var[i]]]))
            formals(mytr) <- eval(myalist)
            transform(x,var[i],post=FALSE) <- mytr

        }
    }
    x$attributes$liability[var] <- liability
    x$attributes$ordinal[var] <- TRUE
    if (!missing(labels)) {
        if (length(var)==1) labels <- list(labels)
        x$attributes$labels[var] <- labels
    }
    x$attributes$nordinal[var] <- K
    x$attributes$normal[var] <- FALSE
    covfix(x,var,NULL,exo=exo) <- 1
    if (is.null(index(x))) index(x) <- reindex(x)
    return(x)
}
