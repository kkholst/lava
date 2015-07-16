##' Wrapper function for mclapply
##'
##' @export
##' @param x function or 'sim' object
##' @param R Number of replications
##' @param f Optional function (i.e., if x is a matrix)
##' @param colnames Optional column names
##' @param messages Messages
##' @param mc.cores Number of cores to use
##' @param blocksize Split computations in blocks
##' @param ... Additional arguments to mclapply
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
##' 
##' val <- sim(onerun,R=10,b0=1,messages=0,mc.cores=1)
##' val
##' val <- sim(val,R=40,b0=1,mc.cores=1) ## append results
##' 
##' summary(val,estimate=c(1,1),confint=c(3,4,6,7),true=c(1,1))
##' summary(val,estimate=c(1,1),se=c(2,5),names=c("Model","Sandwich"))
##' 
##' if (interactive()) {
##'     plot(val,estimate=1,c(2,5),true=1,names=c("Model","Sandwich"))
##'     plot(val,estimate=c(1,1),c(2,5),true=c(1,1),names=c("Model","Sandwich"))
##' }
sim.default <- function(x=NULL,R=100,f=NULL,colnames=NULL,messages=1L,mc.cores=parallel::detectCores(),blocksize=2L*mc.cores,...) {
    requireNamespace("parallel",quietly=TRUE)
    olddata <- NULL
    if (inherits(x,c("data.frame","matrix"))) olddata <- x
    if (inherits(x,"sim")) {
        x <- attr(x,"f")
        if (!is.null(f)) x <- f
    } else {
        if (!is.null(f)) x <- f
        if (!is.function(x)) stop("Expected a function or 'sim' object.")
    }
    if (is.null(x)) stop("Must give new function argument 'f'.")
    res <- val <- NULL
    mycall <- match.call()
    on.exit({
        if (messages>0) close(pb)
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
        return(res)
    })
    nfolds <- max(1,round(R/blocksize))
    idx <- split(1:R,sort((1:R)%%nfolds))
    idx.done <- 0
    count <- 0
    if (messages>0) pb <- txtProgressBar(style=3,width=40)
    for (ii in idx) {
        count <- count+1
        val <- parallel::mclapply(ii,x,mc.cores=mc.cores,...)
        if (messages>0) ##getTxtProgressBar(pb)<(i/R)) {
            setTxtProgressBar(pb, count/length(idx))

        if (is.null(res)) {
            res <- matrix(NA,ncol=length(val[[1]]),nrow=R)
        }
        res[ii,] <- Reduce(rbind,val)
        idx.done <- max(ii)
    }
}

##' @export
"[.sim" <- function (x, i, j, drop = FALSE) {
    atr <- attributes(x)
    class(x) <- "matrix"
    x <- NextMethod("[",drop=FALSE)
    atr.keep <- "call"
    if (missing(j)) atr.keep <- c(atr.keep,"f")
    attributes(x)[atr.keep] <- atr[atr.keep]
    class(x) <- c("sim","matrix")
    x
}

##' @export
print.sim <- function(x,...) {
    attr(x,"f") <- attr(x,"call") <- NULL
    class(x) <- "matrix"
    print(x,...)
}

##' @export
plot.sim <- function(x,estimate=NULL,se=NULL,true=NULL,
                     names=NULL,
                     auto.layout=TRUE,
                     type="p",
                     ask=grDevices::dev.interactive(),                     
                     col=rgb(.5,.5,.5),pch=16,cex=0.5,
                     legend=colnames(x),
                     density.ylim,
                     density.plot=TRUE,running.mean=TRUE,
                     scatter.plot=TRUE,
                     density.lty=2:9,
                     density.col=rgb(.5,.5,.5),
                     xlab="",...) {

    if (is.null(estimate)) {
        av <- apply(x[,drop=FALSE],2,function(z) cumsum(z)/seq(length(z)))
        matplot(av,type="l",main="Running Average",...)
        if (!is.null(true)) abline(h=true,lty=2,...)
        if (!is.null(legend))
            graphics::legend("bottom",legend=legend,bg="white",...)
    }
    
    K <- length(estimate)
    est <- tru <- c()
    if (length(se)>0) {
        if (K==1 && !is.list(se))
            se <- list(se)
        else se <- as.list(se)
    } else {
        est <- estimate; tru <- true
    }
    for (i in seq_along(estimate)) {
        est <- c(est,list(rep(estimate[i],length(se[[i]]))))
        if (!is.null(true)) tru <- c(tru,list(rep(true[i],length(se[[i]]))))
    }
    ss <- summary(x,estimate=unlist(est),se=unlist(se),true=unlist(tru),names=names)
    oldpar <- NULL
    on.exit({
        par(oldpar)
        return(ss)
    })

    if (auto.layout) {
        nc <- (scatter.plot || running.mean) + density.plot
        nr <- min(6,K)
        oma.multi = c(3, 0, 3, 0)
        mar.multi = c(1, 4.1, 1, 1)
        oldpar <- par(mar=mar.multi, oma=oma.multi,
                      mfrow=c(nr,nc),
                      cex.axis=0.5,las=1,
                      ask=FALSE)
    }
    
    if (!missing(density.ylim)) density.ylim <- rep(density.ylim,length.out=K)
    for (i in seq(K)) {
        ii <- estimate[i]
        y <- as.vector(x[,ii])
        if (scatter.plot || running.mean) 
            graphics::plot(y,ylab=colnames(x)[ii],col=col,cex=cex,pch=pch,type=type)
        if (running.mean) {
            lines(cumsum(y)/seq_along(y))
            if (!is.null(true))
                abline(h=true[i],lty=density.lty[1])
        }
        if (density.plot) {
            dy <- stats::density(y)
            if (missing(density.ylim)) {
                density.ylim0 <- c(0,max(dy$y)*1.5)
            } else {
                density.ylim0 <- density.ylim[j]
            }
            graphics::plot(dy,main="",xlab=xlab,ylim=density.ylim0); graphics::rug(y,col=col)
        }
        if (!is.null(se)) {
            se.pos <- match(se[[i]],unlist(se))
            se.col <- rep(density.col,length.out=length(se.pos))
            se.lty <- rep(density.lty,length.out=length(se.pos))
            for (j in seq_along(se.pos)) {
                graphics::curve(dnorm(x,mean=ss["Mean",i],sd=ss["SE",se.pos[j]]),lty=se.lty[j],col=se.col[j],add=TRUE)
            }
            graphics::legend("topright",colnames(ss)[se.pos],col=se.col,lty=se.lty,cex=0.6)
        }
        if (i==1 && ask) par(ask=ask)
    }
}

##' @export
summary.sim <- function(object,estimate=NULL,se=NULL,confint=NULL,true=NULL,fun,names=NULL,...) {
    if (missing(fun)) fun <- function(x) {
        pp <- c(.025,.5,.975)
        res <- c(mean(x,na.rm=TRUE),sd(x,na.rm=TRUE),quantile(x,c(0,pp,1),na.rm=TRUE),
                 mean(is.na(x)))
        names(res) <- c("Mean","SD","Min",paste0(pp*100,"%"),"Max","Missing")
        res
    }
    if (is.null(estimate) && is.null(confint)) return(apply(object,2,fun))

    est <- apply(object[,estimate,drop=FALSE],2,
                 function(x) c(Mean=mean(x,na.rm=TRUE),Missing=mean(is.na(x)),SD=sd(x,na.rm=TRUE)))
    if (!is.null(true)) {
        if (length(true)!=length(estimate)) stop("'true' should be of same length as 'estimate'.")
        est <- rbind(rbind(True=true),rbind(Bias=true-est["Mean",]),
                     rbind(RMSE=((true-est["Mean",])^2+(est["SD",])^2)^.5),
                     est)
    }
    if (!is.null(se)) {
        if (length(se)!=length(estimate)) stop("'se' should be of same length as 'estimate'.")
        est <- rbind(est, SE=apply(object[,se,drop=FALSE],2,
                                  function(x) c(mean(x,na.rm=TRUE))))
        est <- rbind(est,"SE/SD"=est["SE",]/est["SD",])

    }
    if (!is.null(confint)) {
        if (length(confint)!=2*length(estimate)) stop("'confint' should be of length 2*length(estimate).")
        Coverage <- c()
        for (i in seq_along(estimate)) {
            Coverage <- c(Coverage,
                          mean((object[,confint[2*(i-1)+1]]<true[i]) & (object[,confint[2*i]]>true[i]),na.rm=TRUE))
        }
        est <- rbind(est,Coverage=Coverage)
    }
    if (!is.null(names)) colnames(est) <- names
    return(est)
}
