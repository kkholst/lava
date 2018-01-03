##' Wrapper function for mclapply
##'
##' @export
##' @param x function or 'sim' object
##' @param R Number of replications or data.frame with parameters
##' @param f Optional function (i.e., if x is a matrix)
##' @param colnames Optional column names
##' @param messages Messages
##' @param mc.cores Number of cores to use
##' @param cl (optional) cluster to use for parallelization
##' @param blocksize Split computations in blocks
##' @param type type=0 is an alias for messages=1,mc.cores=1,blocksize=R
##' @param seed (optional) Seed (needed with cl=TRUE)
##' @param args (optional) list of named arguments passed to (mc)mapply
##' @param iter If TRUE the iteration number is passed as first argument to (mc)mapply
##' @param ... Additional arguments to (mc)mapply
##' @aliases sim.default
##' @seealso summary.sim plot.sim print.sim
##' @examples
##' m <- lvm(y~x+e)
##' distribution(m,~y) <- 0
##' distribution(m,~x) <- uniform.lvm(a=-1.1,b=1.1)
##' transform(m,e~x) <- function(x) (1*x^4)*rnorm(length(x),sd=1)
##'
##' onerun <- function(iter=NULL,...,n=2e3,b0=1,idx=2) {
##'     d <- sim(m,n,p=c("y~x"=b0))
##'     l <- lm(y~x,d)
##'     res <- c(coef(summary(l))[idx,1:2],
##'              confint(l)[idx,],
##'              estimate(l,only.coef=TRUE)[idx,2:4])
##'     names(res) <- c("Estimate","Model.se","Model.lo","Model.hi",
##'                     "Sandwich.se","Sandwich.lo","Sandwich.hi")
##'     res
##' }
##' val <- sim(onerun,R=10,b0=1,messages=0,mc.cores=1)
##' val
##'
##' val <- sim(val,R=40,b0=1,mc.cores=1) ## append results
##' summary(val,estimate=c(1,1),confint=c(3,4,6,7),true=c(1,1))
##'
##' summary(val,estimate=c(1,1),se=c(2,5),names=c("Model","Sandwich"))
##' summary(val,estimate=c(1,1),se=c(2,5),true=c(1,1),names=c("Model","Sandwich"),confint=TRUE)
##'
##' if (interactive()) {
##'     plot(val,estimate=1,c(2,5),true=1,names=c("Model","Sandwich"),polygon=FALSE)
##'     plot(val,estimate=c(1,1),se=c(2,5),main=NULL,
##'          true=c(1,1),names=c("Model","Sandwich"),
##'          line.lwd=1,density.col=c("gray20","gray60"),
##'          rug=FALSE)
##'     plot(val,estimate=c(1,1),se=c(2,5),true=c(1,1),
##'          names=c("Model","Sandwich"))
##' }
##'
##' f <- function(a=1,b=1) {
##'   rep(a*b,5)
##' }
##' R <- Expand(a=1:3,b=1:3)
##' sim(f,R,type=0)
##' sim(function(a,b) f(a,b), 3, args=c(a=5,b=5),type=0)
##' sim(function(iter=1,a=5,b=5) iter*f(a,b), type=0, iter=TRUE, R=5)
sim.default <- function(x=NULL,R=100,f=NULL,colnames=NULL,
                messages=lava.options()$messages,
                mc.cores,blocksize=2L*mc.cores,
                cl,type=1L,seed=NULL,args=list(),iter=FALSE,...) {
    stm <- proc.time()
    oldtm <- rep(0,5)
    if (missing(mc.cores) || .Platform$OS.type=="windows") {
        if (.Platform$OS.type=="windows") { ## Disable parallel processing on windows
            mc.cores <- 1L
        } else {
            mc.cores <- getOption("mc.cores",parallel::detectCores())
        }
    }
    if (type==0L) {
        mc.cores <- 1L
        if (inherits(R,c("matrix","data.frame")) || length(R)>1) {
            blocksize <- NROW(R)
        } else {
            blocksize <- R
        }
        messages <- 0
    }
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    if (mc.cores>1L || !missing(cl)) requireNamespace("parallel",quietly=TRUE)
    newcl <- FALSE
    if (!missing(cl) && is.logical(cl)) {
        if (.Platform$OS.type=="windows" || TRUE) { ## Don't fork processes on windows
            cl <- NULL
            mc.cores <- 1
        } else {
            if (cl) {
                cl <- parallel::makeForkCluster(mc.cores)
                if (!is.null(seed)) parallel::clusterSetRNGStream(cl,seed)
                newcl <- TRUE
            }
        }
    }
    olddata <- NULL
    dots <- list(...)
    mycall <- match.call(expand.dots=FALSE)
    if (inherits(x,c("data.frame","matrix"))) olddata <- x
    if (inherits(x,"sim")) {
        oldtm <- attr(x,"time")
        oldcall <- attr(x,"call")
        x <- attr(x,"f")
        if (!is.null(f)) x <- f
        ex <- oldcall[["..."]]
        for (nn in setdiff(names(ex),names(dots))) {
            dots[[nn]] <- ex[[nn]]
            val <- list(ex[[nn]]); names(val) <- nn
            mycall[["..."]] <- c(mycall[["..."]],list(val))
        }

    } else {
        if (!is.null(f)) x <- f
        if (!is.function(x)) stop("Expected a function or 'sim' object.")
    }
    if (is.null(x)) stop("Must give new function argument 'f'.")
    res <- val <- NULL
    on.exit({
        if (messages>0) close(pb)
        if (newcl) parallel::stopCluster(cl)
        if (is.null(colnames) && !is.null(val)) {
            if (is.matrix(val[[1]])) {
                colnames <- base::colnames(val[[1]])
            } else {
                colnames <- names(val[[1]])
            }
        }
        base::colnames(res) <- colnames
        if (!is.null(olddata)) res <- rbind(olddata,res)
        attr(res,"call") <- mycall
        attr(res,"f") <- x
        class(res) <- c("sim","matrix")
        if (idx.done<R) {
            res <- res[seq(idx.done),,drop=FALSE]
        }
        attr(res,"time") <- proc.time()-stm+oldtm

        return(res)
    })
    parval_provided <- FALSE
    if (inherits(R,c("matrix","data.frame")) || length(R)>1) {
        parval_provided <- TRUE
        parval <- as.data.frame(R)
        if (is.vector(R)) names(parval) <- NULL
        else if (inherits(R,c("matrix","data.frame"))) names(parval) <- colnames(R)
        R <- NROW(parval)
    } else {
        parval <- as.data.frame(1:R)
        names(parval) <- NULL
    }
    nfolds <- max(1,round(R/blocksize))
    idx <- split(1:R,sort((1:R)%%nfolds))
    idx.done <- 0
    count <- 0
    if (messages>0) pb <- txtProgressBar(style=lava.options()$progressbarstyle,width=40)
    time <- c()
    robx <- function(iter__,...) tryCatch(x(...),error=function(e) NA)
    if (iter) formals(robx)[[1]] <- NULL
    for (ii in idx) {
        count <- count+1
        if (!missing(cl) && !is.null(cl)) {            
            pp <- c(as.list(parval[ii,,drop=FALSE]),dots,list(cl=cl,fun=robx,SIMPLIFY=FALSE),args)
        } else {
            pp <- c(as.list(parval[ii,,drop=FALSE]),dots,list(mc.cores=mc.cores,FUN=robx,SIMPLIFY=FALSE),args)
        }
        ##if (!iter & !parval_provided) pp[[1]] <- NULL
        if (mc.cores>1) {
            if (!missing(cl) && !is.null(cl)) {
                val <- do.call(parallel::clusterMap,pp)
            } else {
                val <- do.call(parallel::mcmapply,pp)
            }
        } else {
            pp$mc.cores <- NULL
            val <- do.call(mapply,pp)
        }
        if (messages>0)
            setTxtProgressBar(pb, count/length(idx))
        if (is.null(res)) {
            ##res <- array(NA,dim=c(R,dim(val[[1]])),dimnames=c(list(NULL),dimnames(val[[1]]),NULL))
            res <- matrix(NA,ncol=length(val[[1]]),nrow=R)
        }
        res[ii,] <- Reduce(rbind,val)
        ##rr <- abind::abind(val,along=length(dim(res)))
        ##res[ii,] <- abind(val,along=length(dim(res)))
        idx.done <- max(ii)
    }
}

##' @export
"[.sim" <- function (x, i, j, drop = FALSE) {
    atr <- attributes(x)
    if (!is.null(dim(x))) {
        class(x) <- "matrix"
    } else {
        class(x) <- class(x)[-1]
    }
    x <- NextMethod("[",drop=drop)
    atr.keep <- c("call","time")
    if (missing(j)) atr.keep <- c(atr.keep,"f")
    attributes(x)[atr.keep] <- atr[atr.keep]
    if (!drop) class(x) <- c("sim",class(x))
    x
}


Time <- function(sec,print=FALSE,...) {
    h <- sec%/%3600
    m0 <- (sec%%3600)
    m <- m0%/%60
    s <- m0%%60
    res <- c(h=h,m=m,s=s)
    if (print) {
        if (h>0) cat(h,"h ",sep="")
        if (m>0) cat(m,"m ",sep="")
        cat(s,"s",sep="")
        return(invisible(res))
    }
    return(res)
}

Print <- function(x,n=5,digits=max(3,getOption("digits")-3),...) {
    mat <- !is.null(dim(x))
    if (!mat) {
        x <- cbind(x)
        colnames(x) <- ""
    }
    if (is.null(rownames(x))) {
        rownames(x) <- seq(nrow(x))
    }
    sep <- rbind("---"=rep('',ncol(x)))
    if (n<1) {
        print(x,quote=FALSE,digits=digits,...)
    } else {
        ## hd <- base::as.matrix(base::format(utils::head(x,n),digits=digits,...))
        ## tl <- base::as.matrix(base::format(utils::tail(x,n),digits=digits,...))
        ## print(rbind(hd,sep,tl),quote=FALSE,...)
        if (NROW(x)<=(2*n)) {
            hd <- base::format(utils::head(x,2*n),digits=digits,...)
            print(hd, quote=FALSE,...)
        } else {
            hd <- base::format(utils::head(x,n),digits=digits,...)
            tl <- base::format(utils::tail(x,n),digits=digits,...)
            print(rbind(base::as.matrix(hd),sep,base::as.matrix(tl)),
                  quote=FALSE,...)
        }
    }
    invisible(x)
}

##' @export
print.sim <- function(x,...) {
    s <- summary(x,minimal=TRUE,...)
    attr(x,"f") <- attr(x,"call") <- NULL
    if (!is.null(dim(x))) {
        class(x) <- "matrix"
    }
    Print(x,...)
    cat("\n")
    print(s,extra=FALSE,...)
    return(invisible(x))
}



##' @export
print.summary.sim <- function(x,group=list(c("^mean$","^sd$","^se$","^se/sd$","^coverage"),
                                   c("^min$","^[0-9.]+%$","^max$"),
                                   c("^na$","^missing$"),
                                   c("^true$","^bias$","^rmse$")),
                      lower.case=TRUE,
                      na.print="",
                      digits = max(3, getOption("digits") - 2),
                      quote=FALSE,
                      time=TRUE,
                      extra=TRUE,
                      ...) {
    if (extra) {
        cat(attr(x,"n")," replications",sep="")
        if (time && !is.null(attr(x,"time"))) {
            cat("\t\t\t\t\tTime: ")
            Time(attr(x,"time")["elapsed"],print=TRUE)
        }
        cat("\n\n")
    }

    nn <- rownames(x)
    if (lower.case)  nn <- tolower(nn)
    gg <- lapply(group,
                 function(x) unlist(lapply(x,function(v) grep(v,nn))))
    gg <- c(gg,list(setdiff(seq_along(nn),unlist(gg))))

    x0 <- c()
    ng <- length(gg)
    for (i in seq(ng)) {
        x0 <- rbind(x0, x[gg[[i]],,drop=FALSE],
        { if(i<ng && length(gg[[i+1]])>0) NA})
    }

    print(structure(x0,class="matrix")[,,drop=FALSE],digits=digits,quote=quote,na.print=na.print,...)
    if (extra) cat("\n")
    invisible(x)
}


##' Summary method for 'sim' objects
##'
##' Summary method for 'sim' objects
##' @export
##' @export summary.sim
##' @param object sim object
##' @param estimate (optional) columns with estimates
##' @param se (optional) columns with standard error estimates
##' @param confint (optional) list of pairs of columns with confidence limits
##' @param true (optional) vector of true parameter values
##' @param fun (optional) summary function
##' @param names (optional) names of 
##' @param unique.names if TRUE, unique.names will be applied to column names
##' @param minimal if TRUE, minimal summary will be returned
##' @param level confidence level 
##' @param quantiles quantiles
##' @param ... additional levels to lower-level functions
summary.sim <- function(object,estimate=NULL,se=NULL,
                confint=!is.null(se)&&!is.null(true),true=NULL,
                fun,names=NULL,unique.names=TRUE,minimal=FALSE,
                level=0.95,quantiles=c(.025,0.5,.975),...) {
    if (minimal) {
        fun <- function(x,se,confint,...) {
            res <- c(Mean=mean(x,na.rm=TRUE),
                    SD=sd(x,na.rm=TRUE))
            if (!missing(se) && !is.null(se)) {
                res <- c(res, c(SE=mean(se,na.rm=TRUE)))
                res <- c(res, c("SE/SD"=res[["SE"]]/res[["SD"]]))
            }            
            return(res)
        }
    }
    mfun <- function(x,...) {
        res <- c(mean(x,na.rm=TRUE),
                 sd(x,na.rm=TRUE),
                 quantile(x,c(0,quantiles,1),na.rm=TRUE),
                 mean(is.na(x)))
        names(res) <- c("Mean","SD","Min",paste0(quantiles*100,"%"),"Max","Missing")
        res
    }
    tm <- attr(object,"time")
    N <- max(length(estimate),length(se),length(true))
    if (!is.null(estimate)) estimate <- rep(estimate,length.out=N)
    if (!is.null(se)) se <- rep(se,length.out=N)
    if (!is.null(true)) true <- rep(true,length.out=N)
    
    if (!is.null(estimate) && is.character(estimate)) {
        estimate <- match(estimate,colnames(object))
    }
    if (!missing(fun)) {
        if (!is.null(estimate)) m.est <- object[,estimate,drop=FALSE]
        else m.est <- object
        m.se <- NULL
        if (!is.null(se)) m.se <- object[,se,drop=FALSE]
        m.ci <- NULL
        if (!is.null(confint)) m.ci <- object[,confint,drop=FALSE]
        res <- lapply(seq(ncol(m.est)),
                      function(i,...) fun(m.est[,i,drop=TRUE],se=m.se[,i,drop=TRUE],confint=m.ci[,1:2+(i-1)*2],...,INDEX=i),...)
        res <- matrix(unlist(res),nrow=length(res[[1]]),byrow=FALSE)
        if (is.null(dim(res))) {
            res <- rbind(res)
        }
        if (is.null(rownames(res))) {
            rownames(res) <- names(fun(0,m.se,m.ci,INDEX=1,...))
            if (is.null(rownames(res))) rownames(res) <- rep("",nrow(res))
        }
        if (is.null(colnames(res))) {
            colnames(res) <- colnames(m.est)
        }
        return(structure(res,
                    n=NROW(object),
                    time=tm,
                    class=c("summary.sim","matrix")))
    }

    if (!is.null(estimate)) {
        est <- apply(object[,estimate,drop=FALSE],2,mfun)
    } else {
        est <- apply(object,2,mfun)
    }

    if (!is.null(true)) {
        if (length(true)!=length(estimate)) {
            ##stop("'true' should be of same length as 'estimate'.")
            true <- rep(true,length.out=length(estimate))
        }
        est <- rbind(est,
                     rbind(True=true),rbind(Bias=est["Mean",]-true),
                     rbind(RMSE=((est["Mean",]-true)^2+(est["SD",])^2)^.5)
                     )
    }
    if (!is.null(se)) {
        if (is.character(se)) {
            se <- match(se,colnames(object))
        }
        if (length(se)!=length(estimate)) stop("'se' should be of same length as 'estimate'.")
        est <- rbind(est, SE=apply(object[,se,drop=FALSE],2,
                                  function(x) c(mean(x,na.rm=TRUE))))
        est <- rbind(est,"SE/SD"=est["SE",]/est["SD",])

    }
    if (!is.null(confint) && (length(confint)>1 || confint)) {
        if (is.character(confint)) {
            confint <- match(confint,colnames(object))
        }
        if (length(confint)==1 && confint) {
            if (is.null(se)) stop("Supply confidence limits or SE")
            confint <- c()
            pos <- ncol(object)
            for (i in seq_along(estimate)) {
                z <- 1-(1-level)/2
                CI <- cbind(object[,estimate[i]]-qnorm(z)*object[,se[i]],
                            object[,estimate[i]]+qnorm(z)*object[,se[i]])
                colnames(CI) <- NULL
                object <- cbind(object,CI)
                confint <- c(confint,pos+1:2)
                pos <- pos+2
            }
        }
        if (length(confint)!=2*length(estimate)) stop("'confint' should be of length 2*length(estimate).")
        Coverage <- c()
        for (i in seq_along(estimate)) {
            Coverage <- c(Coverage,
                          mean((object[,confint[2*(i-1)+1]]<true[i]) & (object[,confint[2*i]]>true[i]),na.rm=TRUE))
        }
        est <- rbind(est,Coverage=Coverage)
    }
    if (!is.null(names)) {
         if (length(names)<ncol(est)) {
            uest <- unique(estimate)
            names <- names[match(estimate,uest)]
        }
        colnames(est) <- names

    }
    if (unique.names && !is.null(colnames(est))) {
        colnames(est) <- make.unique(colnames(est))
    }

    return(structure(est,
                     n=NROW(object),
                     time=tm,
                     class=c("summary.sim","matrix")))
}
