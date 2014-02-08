##' Estimation of functional of parameters 
##'
##' Estimation of functional of parameters.
##' Wald tests, robust standard errors, cluster robust standard errors,
##' LRT (when \code{f} is not a function)...
##' @param x model object (\code{glm}, \code{lvmfit}, ...)
##' @param f transformation of model parameters and (optionally) data, or contrast matrix (or vector)
##' @param data \code{data.frame}
##' @param id (optional) id-variable corresponding to iid decomposition of model parameters. 
##' @param stack 
##' @param subset (optional) subset of data.frame on which to condition (logical expression or variable name)
##' @param score.deriv (optional) derivative of mean score function
##' @param level Level of confidence limits
##' @param iid If TRUE the iid decompositions are also returned (extract with \code{iid} method)
##' @param contrast (optional) Contrast matrix for final Wald test
##' @param null (optional) Null hypothesis to test 
##' @param vcov (optional) covariance matrix of parameter estimates (e.g. Wald-test)
##' @param coef (optional) parameter coefficient
##' @param print (optional) print function
##' @param ... additional arguments to lower level functions
##' @export
##' @examples
##' 
##' ## Simulation from logistic regression model
##' m <- lvm(y~x+z);
##' distribution(m,y~x) <- binomial.lvm("logit")
##' d <- sim(m,1000)
##' g <- glm(y~z+x,data=d,family=binomial())
##' g0 <- glm(y~1,data=d,family=binomial())
##' 
##' ## LRT
##' estimate(g,g0)
##' 
##' ## Plain estimates (robust standard errors)
##' estimate(g)
##' 
##' ## Testing contrasts
##' estimate(g,null=0)
##' estimate(g,rbind(c(1,1,0),c(1,0,2)))
##' estimate(g,rbind(c(1,1,0),c(1,0,2)),null=c(1,2))
##' estimate(g,2:3) ## same as rbind(c(0,1,0),c(0,0,1))
##' 
##' ## Transformations
##' estimate(g,function(p) p[1]+p[2])
##' 
##' ## Multiple parameters
##' e <- estimate(g,function(p) c(p[1]+p[2],p[1]*p[2]))
##' e
##' vcov(e)
##' 
##' ## Label new parameters
##' estimate(g,function(p) list("a1"=p[1]+p[2],"b1"=p[1]*p[2]))
##' 
##' ## Marginalize
##' f <- function(p,data)
##'   list(p0=lava:::expit(p[1] + p[3]*data[,"z"]),
##'        p1=lava:::expit(p[1] + p[2] + p[3]*data[,"z"]))
##' e <- estimate(g, f)
##' e
##' estimate(e,diff)
##' estimate(e,cbind(1,1))
##'
##' ## Clusters and subset (conditional marginal effects)
##' d$id <- rep(seq(nrow(d)/4),each=4)
##' estimate(g,function(p,data) list(p0=lava:::expit(p[1] + p["z"]*data[,"z"])), subset=z>0, id=d$id)
##' 
##' @method estimate default
##' @S3method estimate default
estimate.default <- function(x,f=NULL,data=model.frame(x),id,stack=TRUE,subset,
                             score.deriv,level=0.95,iid=TRUE,
                             contrast,null,vcov,coef,print=NULL,...) {
    if (!is.null(f) && !is.function(f)) {
        if (!(is.matrix(f) | is.vector(f))) return(compare(x,f,...))
        contrast <- f; f <- NULL
    }
    if (is.matrix(x) || is.vector(x)) contrast <- x
    alpha <- 1-level
    alpha.str <- paste(c(alpha/2,1-alpha/2)*100,"",sep="%")
    nn <- NULL
    if (missing(vcov)) { ## If user supplied vcov, then don't estimate IC
        if (missing(score.deriv)) {
            if (!is.logical(iid)) {
                iidtheta <- iid
                iid <- TRUE
            } else {
                suppressWarnings(iidtheta <- iid(x))
            }
        } else {
            suppressWarnings(iidtheta <- iid(x,score.deriv=score.deriv))
        }
    }  else { iidtheta <- NULL }
    if (!missing(subset)) {
        e <- substitute(subset)
        subset <- eval(e, data, parent.frame())
        if (is.character(subset)) subset <- data[,subset]
        if (is.numeric(subset)) subset <- subset==1
    }
    if (!missing(id)) {
        if (is.null(iidtheta)) stop("'iid' method needed")
        nprev <- nrow(iidtheta)
        e <- substitute(id)        
        id <- eval(e, data, parent.frame())
        if (is.logical(id) && length(id)==1) {
            id <- if(is.null(iidtheta)) seq(nrow(data)) else seq(nprev)
            stack <- FALSE
        }
        if (is.character(id) && length(id)==1) id <- data[,id,drop=TRUE]
        if (!is.null(iidtheta)) {
            if (length(id)!=nprev) stop("Dimensions of i.i.d decomposition and 'id' does not agree")
        } else {
            if (length(id)!=nrow(data)) stop("Dimensions of 'data' and 'id' does not agree")
        }
        if (stack) {
            clidx <- NULL
            if (inherits(try(find.package("mets"),silent=TRUE),"try-error")) {
                iidtheta <- matrix(unlist(by(iidtheta,id,colSums)),byrow=TRUE,ncol=ncol(iidtheta))
            } else {
                clidx <- mets::cluster.index(id)
                iidtheta <- t(rbind(apply(clidx$idclustmat+1,1,function(x) colSums(iidtheta[na.omit(x),,drop=FALSE]))))
            }            
            id <- id[as.vector(clidx$firstclustid)+1]
        }
    } else id <- NULL
    if (!is.null(iidtheta)) rownames(iidtheta) <- id

    if (!is.null(iidtheta) && missing(vcov)) {
        n <- NROW(iidtheta)
        if (is.null(f)) {
            V <- crossprod(iidtheta)
        }
    } else {
        if (!missing(vcov)) {
            V <- vcov
        } else {
            V <- stats::vcov(x)
        }
    }
    if (!missing(coef)) {
        pp <- coef
    } else {
        pp <- stats::coef(x)
    }

    if (!is.null(f)) {
        form <- names(formals(f))
        dots <- ("..."%in%names(form))
        form0 <- setdiff(form,"...")
        parname <- "p"
        parname <- form[1]
        if (length(form0)==1 && !(form0%in%c("object","data"))) {
            ##names(formals(f))[1] <- "p"
            parname <- form0
        }
        if (!is.null(iidtheta)) {
            arglist <- c(list(object=x,data=data,p=pp),list(...))
            names(arglist)[3] <- parname      
        } else {
            arglist <- c(list(object=x,p=pp),list(...))            
            names(arglist)[2] <- parname
        }
        if (!dots) {
            arglist <- arglist[form0]
        }

        newf <- NULL
        if (length(form)==0) {
            arglist <- list(pp)
            newf <- function(p,...) do.call("f",list(p,...))
            newf <- function(...) do.call("f",list(...))
            val <- do.call("f",arglist)
        } else {
            val <- do.call("f",arglist)
            if (is.list(val)) {
                nn <- names(val)
                val <- do.call("cbind",val)
                newf <- function(p,...) do.call("cbind",f(p,...))
                newf <- function(...) do.call("cbind",f(...))
            }
        }
        k <- NCOL(val)
        N <- NROW(val)
        D <- attributes(val)$grad
        if (is.null(D)) {
            D <- numDeriv::jacobian(function(p,...) {
                if (length(form)==0) arglist[[1]] <- p
                else arglist[[parname]] <- p
                if (is.null(newf))
                    return(do.call("f",arglist))
                return(do.call("newf",arglist)) }, pp)      
        }
        if (is.null(iidtheta)) {
            pp <- as.vector(val)
            V <- D%*%V%*%t(D)
        } else {
            if (N<NROW(data) || NROW(data)==0) { ## transformation not depending on data
                pp <- as.vector(val)
                iidtheta <- iidtheta%*%t(D)
                V <- crossprod(iidtheta)
            } else {
                if (k>1) { ## More than one parameter (and depends on data)
                    if (!missing(subset)) { ## Conditional estimate
                        val <- apply(val,2,function(x) x*subset)
                    }
                    D0 <- matrix(nrow=k,ncol=length(pp))
                    for (i in seq_len(k)) {
                        D1 <- D[seq(N)+(i-1)*N,,drop=FALSE]
                        if (!missing(subset)) ## Conditional estimate
                            D1 <- apply(D1,2,function(x) x*subset)
                        D0[i,] <- colMeans(D1)
                    }
                    D <- D0
                    iid2 <- iidtheta%*%t(D)        
                } else { ## Single parameter
                    if (!missing(subset)) { ## Conditional estimate
                        val <- val*subset
                        D <- apply(rbind(D),2,function(x) x*subset)
                    }
                    D <- colMeans(rbind(D))
                    iid2 <- iidtheta%*%D
                }
                pp <- as.vector(colMeans(cbind(val)))
                iid1 <- (cbind(val)-rbind(pp)%x%cbind(rep(1,N)))/N
                if (!missing(id)) {
                    if (is.null(clidx)) 
                        iid1 <- matrix(unlist(by(iid1,id,colSums)),byrow=TRUE,ncol=ncol(iid1))
                    else {
                        iid1 <- t(rbind(apply(clidx$idclustmat+1,1,function(x) colSums(iid1[na.omit(x),,drop=FALSE]))))
                    }
                }
                if (!missing(subset)) { ## Conditional estimate
                    phat <- mean(subset)
                    iid3 <- cbind(-1/phat^2 * (subset-phat)/n)
                    if (!missing(id)) {
                        if (is.null(clidx))
                            iid3 <- matrix(unlist(by(iid3,id,colSums)),byrow=TRUE,ncol=ncol(iid3))
                        else
                            iid3 <- cbind(apply(clidx$idclustmat+1,1,function(x) sum(iid3[na.omit(x)])))
         
                    }
                    iidtheta <- (iid1+iid2)/phat + rbind(pp)%x%iid3
                    pp <- pp/phat
                    V <- crossprod(iidtheta)
                } else {
                    if (nrow(iid1)!=nrow(iid2)) {
                        message("Assuming independence between model iid decomposition and new data frame")
                        V <- crossprod(iid1) + crossprod(iid2)
                    } else {
                        iidtheta <- iid1+iid2
                        V <- crossprod(iidtheta)
                    }
                }
            }            
        }
    }

    if (length(pp)==1) res <- rbind(c(pp,diag(V)^0.5)) else res <- cbind(pp,diag(V)^0.5)
    if (is.null(f) || missing(null))
        res <- cbind(res,res[,1]-qnorm(1-alpha/2)*res[,2],res[,1]+qnorm(1-alpha/2)*res[,2],(pnorm(abs(res[,1]/res[,2]),lower.tail=FALSE)*2))
    else 
        res <- cbind(res,res[,1]-qnorm(1-alpha/2)*res[,2],res[,1]+qnorm(1-alpha/2)*res[,2],(pnorm(abs(res[,1]-null)/res[,2],lower.tail=FALSE))*2)
    colnames(res) <- c("Estimate","Std.Err",alpha.str,"P-value")
    if (!is.null(nn)) {
        rownames(res) <- nn
    } else {
        nn <- attributes(res)$varnames
        if (!is.null(nn)) rownames(res) <- nn
        if (is.null(rownames(res))) rownames(res) <- paste("p",seq(nrow(res)),sep="")
    }
    coefs <- res[,1,drop=TRUE]; names(coefs) <- rownames(res)
    res <- structure(list(coef=coefs,coefmat=res,vcov=V, iid=NULL, print=print, id=id),class="estimate")
    if (iid) res$iid <- iidtheta
    if (is.null(f) && (!missing(contrast) | !missing(null))) {
        p <- length(res$coef)    
        if (missing(contrast)) contrast <- diag(p)
        if (missing(null)) null <- 0
        if (is.vector(contrast)) {
            if (length(contrast)==p) contrast <- rbind(contrast)
            else {
                cont <- contrast
                contrast <- diag(nrow=p)[cont,,drop=FALSE]
            }
        }
        cc <- compare(res,contrast=contrast,null=null,vcov=V,level=level)
        res <- structure(c(res, list(compare=cc)),class="estimate")
        res$coefmat <- with(cc, cbind(estimate,
                                      (pnorm(abs(estimate[,1]-null)/estimate[,2],lower.tail=FALSE)*2)));
        colnames(res$coefmat)[5] <- "P-value"
        rownames(res$coefmat) <- cc$cnames
        res$compare$estimate <- NULL
        res$coef <- res$compare$coef
        res$vcov <- res$compare$vcov
    }
    return(res)  
}

estimate.glm <- function(x,...) {  
    estimate.default(x,...)
}

##' @S3method print estimate
print.estimate <- function(x,digits=3,width=25,...) {
    if (!is.null(x$print)) {
        x$print(x,...)
        return(invisible(x))
    }
    cc <- x$coefmat
    rownames(cc) <- make.unique(unlist(lapply(rownames(cc),
                                               function(x) toString(x,width=width))))
    print(cc,digits=digits,...)
    if (!is.null(x$compare)) print(x$compare)    
}

##' @S3method vcov estimate
vcov.estimate <- function(object,...) {    
    res <- object$vcov
    nn <- names(coef(object))
    dimnames(res) <- list(nn,nn)
    res
}

##' @S3method coef estimate
coef.estimate <- function(object,...) {
    object$coef
}

##' @S3method iid estimate
iid.estimate <- function(x,...) {
    x$iid
}

##' @S3method model.frame estimate
model.frame.estimate <- function(formula,...) {
    NULL
}

##' @S3method merge estimate
merge.estimate <- function(x,y,...,id,indep=FALSE) {
    objects <- list(x,y, ...)
    coefs <- unlist(lapply(objects,coef))
    names(coefs) <- make.unique(names(coefs))   
    if (missing(id)) {
        id <- lapply(objects,function(x) x$id)
    } else {
        if (length(id)!=length(objects)) stop("Same number of id-elements as model objects expected")
    }
    if (any(unlist(lapply(id,is.null)))) stop("Id needed for each model object")
    ##iid <- Reduce("cbind",lapply(objects,iid))    
    ids <- iidall <- c(); count <- 0
    for (z in objects) {
        count <- count+1
        clidx <- NULL
        id0 <- id[[count]]
        if (inherits(try(find.package("mets"),silent=TRUE),"try-error")) {
            iid0 <- matrix(unlist(by(iid(z),id0,colSums)),byrow=TRUE,ncol=length(coef(z)))
        } else {
            clidx <- mets::cluster.index(id0)
            iid0 <- t(rbind(apply(clidx$idclustmat+1,1,function(x) colSums(iid(z)[na.omit(x),,drop=FALSE]))))
        }
        ids <- c(ids, list(id0[as.vector(clidx$firstclustid)+1]))
        iidall <- c(iidall, list(iid0))
    }
    id <- unique(unlist(ids))
    iid0 <- matrix(0,nrow=length(id),ncol=length(coefs))
    colpos <- 0
    for (i in seq(length(objects))) {
        relpos <- seq(length(coef(objects[[i]])))
        iid0[match(ids[[i]],id),relpos+colpos] <- iidall[[i]]
        colpos <- colpos+tail(relpos,1)
    }
    rownames(iid0) <- id##
    estimate.default(NULL,coef=coefs, stack=FALSE, data=NULL, iid=iid0, id=id)
}

##' @S3method %+% estimate
`%+%.estimate` <- function(x,y) merge(x,y)



