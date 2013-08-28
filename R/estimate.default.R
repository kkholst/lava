0##' Estimation of functional of parameters 
##'
##' Estimation of functional of parameters.
##' Wald tests, robust standard errors, cluster robust standard errors,
##' LRT (when \code{f} is not a function)...
##' @param x model object (\code{glm}, \code{lvmfit}, ...)
##' @param f transformation of model parameters and (optionally) data
##' @param data \code{data.frame}
##' @param id (optional) id-variable corresponding to iid decomposition of model parameters
##' @param id2 (optional) id-variable of data.frame
##' @param score.deriv (optional) derivative of mean score function
##' @param level Level of confidence limits
##' @param iid If TRUE the iid decompositions are returned instead of variance estimates
##' @param contrast (optional) Contrast matrix for final Wald test
##' @param null (optional) Null hypothesis to test 
##' @param vcov (optional) covariance matrix of parameter estimates (e.g. Wald-test)
##' @param coef (optional) parameter coefficient
##' @param ... additional arguments to lower level functions
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
##' estimate(g,contrast=rbind(c(1,1,0),c(1,0,2)))
##' estimate(g,contrast=rbind(c(1,1,0),c(1,0,2)),null=c(1,2))
##' estimate(g,contrast=2:3) ## same as rbind(c(0,1,0),c(0,0,1))
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
##' estimate(e,contrast=cbind(1,1))
##' @method estimate default
##' @S3method estimate default
estimate.default <- function(x,f,data=model.frame(x),id,id2,
                             score.deriv,level=0.95,iid=FALSE,
                             contrast,null,vcov,coef,...) {
    if (!missing(f) && !is.function(f)) return(compare(x,f,...))
    alpha <- 1-level
    alpha.str <- paste(c(alpha/2,1-alpha/2)*100,"",sep="%")
    nn <- NULL
    if (missing(score.deriv)) {
        suppressWarnings(iid0 <- iid(x))
    } else {
        suppressWarnings(iid0 <- iid(x,score.deriv=score.deriv))
    }
    if (!missing(id)) {
        if (is.null(iid0)) stop("'iid' method needed")
        if (is.character(id) && length(id)==1) id <- data[,id,drop=TRUE]
        nprev <- nrow(iid0)
        iid0 <- matrix(unlist(by(iid0,id,colSums)),byrow=TRUE,ncol=ncol(iid0))
    }

    if (!is.null(iid0) && missing(vcov)) {
        n <- NROW(iid0)
        if (missing(f)) {
            if (iid) return(iid0)
            V <- crossprod(iid0)
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

    if (!missing(f)) {
        form <- names(formals(f))
        dots <- ("..."%in%names(form))
        form0 <- setdiff(form,"...")
        parname <- "p"
        if (length(form0)==1 && !(form0%in%c("object","data"))) {
            ##names(formals(f))[1] <- "p"
            parname <- form0
        }
        if (!is.null(iid0)) {
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
            newf <- function(p) do.call("f",list(p))
            val <- do.call("f",arglist)
        } else {
            val <- do.call("f",arglist)
            if (is.list(val)) {
                nn <- names(val)
                val <- do.call("cbind",val)
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

        if (is.null(iid0)) {
            pp <- as.vector(val)
            V <- D%*%V%*%t(D)
        } else {
            if (NROW(val)<NROW(data)) { ## Apparently transformation not depending on data
                pp <- as.vector(val)
                iid2 <- iid0%*%t(D)
                V <- crossprod(iid2)
            } else {      
                if (k>1) { ## More than one parameter (and depends on data)
                    D0 <- matrix(nrow=k,ncol=length(pp))
                    for (i in seq_len(k)) {
                        D0[i,] <- colMeans(D[seq(N)+(i-1)*N,,drop=FALSE])
                    }
                    D <- D0
                    iid2 <- iid0%*%t(D)        
                } else {
                    D <- colMeans(rbind(D))
                    iid2 <- iid0%*%D
                }
                pp <- as.vector(colMeans(cbind(val)))
                iid1 <- t(rbind(apply(cbind(val),1,function(x) x-pp)))/length(val)
                if (iid) return(list(iid1,iid2))
                if (N!=n) {
                    ## TODO....
                    ##                    browser()
                    if (missing(id2)) {                        
                        message("Assuming independence between model iid decomposition and new data frame")
                        V <- crossprod(iid1) + crossprod(iid2)
                    } else {
                        if (missing(id)) {
                            ii <- iid2
                            ii[id2,] <- ii[id2,] + iid1
                            V <- crossprod(ii)
                        }
                        else stop("No implementation of id,id2 combination...")

                    }
                } else {
                    V <- crossprod((iid1+iid2))
                }
            }
        }
    }
    
    res <- cbind(pp,diag(V)^0.5)
    if (missing(f) || missing(null))
        res <- cbind(res,res[,1]-qnorm(1-alpha/2)*res[,2],res[,1]+qnorm(1-alpha/2)*res[,2],(1-pnorm(abs(res[,1])/res[,2]))*2)
    else 
        res <- cbind(res,res[,1]-qnorm(1-alpha/2)*res[,2],res[,1]+qnorm(1-alpha/2)*res[,2],(1-pnorm(abs(res[,1]-null)/res[,2]))*2)  
    colnames(res) <- c("Estimate","Std.Err",alpha.str,"P-value")
    if (!is.null(nn)) {
        rownames(res) <- nn
    } else {
        nn <- attributes(res)$varnames
        if (!is.null(nn)) rownames(res) <- nn
        if (is.null(rownames(res))) rownames(res) <- paste("p",seq(nrow(res)),sep="")
    }
    res <- structure(list(coef=res[,1],coefmat=res,vcov=V),class="estimate")
    if (missing(f) && (!missing(contrast) | !missing(null))) {
        p <- length(res$coef)    
        if (missing(contrast)) contrast <- diag(p)
        if (missing(null)) null <- 0
        if (is.vector(contrast)) {
            cont <- contrast
            contrast <- diag(nrow=p)[cont,,drop=FALSE]
        }
        cc <- compare(res,contrast=contrast,null=null,vcov=V)
        res <- structure(c(res, list(compare=cc)),class="estimate")
    }
    return(res)  
}

estimate.glm <- function(x,...) {  
    estimate.default(x,...,var)
}


##' @S3method print estimate
print.estimate <- function(x,digits=3,...) {
    ##  cat("\n")
    print(x$coefmat,digits=digits,...)
    if (!is.null(x$compare)) print(x$compare)    
}

##' @S3method vcov estimate
vcov.estimate <- function(object,...) {
    object$vcov
}

##' @S3method coef estimate
coef.estimate <- function(object,...) {
    object$coef
}


