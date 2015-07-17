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
##'     plot(val,estimate=c(1,1),se=c(2,5),true=c(1,1),
##'          names=c("Model","Sandwich"))
##'     plot(val,estimate=c(1,1),se=c(2,5),main=NULL,
##'          true=c(1,1),names=c("Model","Sandwich"),
##'          polygon=TRUE,line.lwd=1.5,density.col=c("gray20","gray60"),
##'          rug=TRUE)
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
                     byrow=FALSE,
                     type="p",
                     ask=grDevices::dev.interactive(),
                     line.col=1, line.lwd=1,
                     col=rgb(.5,.5,.5),pch=16,cex=0.5,lty=1,
                     legend=colnames(x),
                     legendpos="topleft",
                     polygon=FALSE,
                     rug=!polygon,
                     main,
                     ylab="Estimate",
                     density.ylim,
                     density.xlim,
                     density.plot=TRUE,
                     scatter.plot=TRUE,
                     running.mean=scatter.plot,
                     alpha=0.25,
                     border=NULL,
                     density.lty=1:9,
                     density.col=c("gray20"),
                     density.lwd=1,
                     xlab="",...) {

    if (is.null(estimate)) {
        av <- apply(x[,drop=FALSE],2,function(z) cumsum(z)/seq(length(z)))
        matplot(x,type="p",pch=pch, cex=cex, col=col,...)
        matlines(av,type="l",col=col,lty=lty,...)
        if (!is.null(true)) abline(h=true,lty=2,...)
        if (!is.null(legend))
            graphics::legend(legendpos,legend=legend,bg="white",col=col,lty=lty,pch=pch,...)
        return(invisible(NULL))
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
        return(invisible(ss))
    })

    if (auto.layout) {
        nc <- (scatter.plot || running.mean) + density.plot
        nr <- min(6,K)
        oma.multi = c(2, 0, 2, 0)
        mar.multi = c(1.5, 4.1, 1, 1)
        oldpar <- par(mar=mar.multi, oma=oma.multi,
                      cex.axis=0.5,las=1,
                      ask=FALSE)
        if (byrow) {
            par(mfrow=c(nr,nc))
        } else {
            par(mfcol=c(nc,nr))
        }
    }
    
    if (!missing(density.ylim)) density.ylim <- rep(density.ylim,length.out=K)
    if (!missing(density.xlim)) density.xlim <- rep(density.xlim,length.out=K)
    if (missing(main)) {
        main <- colnames(ss)
    }
    if (!is.null(main)) main <- rep(main,length.out=K)
    for (i in seq(K)) {
        ii <- estimate[i]
        y <- as.vector(x[,ii])
        if (scatter.plot || running.mean)            
            graphics::plot(y,ylab=ylab,col=col,cex=cex,pch=pch,type=type)
        if (!is.null(main) && !byrow) {
            title(main[i])
        }
        if (running.mean) {
            lines(cumsum(y)/seq_along(y),col=line.col,lwd=line.lwd,lty=lty)
            if (!is.null(true))
                abline(h=true[i],lty=density.lty[1],col=line.col)
        }
        if (density.plot) {
            dy <- stats::density(y)
            if (missing(density.ylim)) {
                density.ylim0 <- c(0,max(dy$y)*1.5)
            } else {
                density.ylim0 <- density.ylim[j]
            }
            if (missing(density.xlim)) {
                density.xlim0 <- range(dy$x)
            } else {
                density.xlim0 <- density.xlim[j]
            }
            graphics::plot(0,0,type="n",main="",ylab="Density",xlab=xlab,ylim=density.ylim0,xlim=density.xlim0)
            if (polygon) {
                with(dy, graphics::polygon(c(x,rev(x)),c(y,rep(0,length(y))),col=Col(density.col[1],alpha=alpha),border=border))
            } else {
                graphics::lines(dy,main="",lty=density.lty[1],col=density.col[1],lwd=density.lwd[1]);
            }
            if (rug) graphics::rug(y,col=col)
            if (!is.null(main) && !(running.mean || scatter.plot)) {
                title(main[i])
            }           
        
            if (!is.null(se)) {
                se.pos <- match(se[[i]],unlist(se))
                se.col <- rep(density.col,length.out=length(se.pos)+1)[-1]
                se.lty <- rep(density.lty,length.out=length(se.pos)+1)[-1]
                se.lwd <- rep(density.lwd,length.out=length(se.pos)+1)[-1]
                xx <- dy$x
                 for (j in seq_along(se.pos)) {
                    if (polygon) {
                        yy <- dnorm(xx,mean=ss["Mean",i],sd=ss["SE",se.pos[j]])
                        graphics::polygon(c(xx,rev(xx)),c(yy,rep(0,length(yy))),col=Col(se.col[1],alpha=alpha),border=border)
                    } else {
                        graphics::curve(dnorm(x,mean=ss["Mean",i],sd=ss["SE",se.pos[j]]),lwd=se.lwd[j],lty=se.lty[j],col=se.col[j],add=TRUE)
                    }
                }
                if (!is.null(legend)) {
                    if (polygon) {
                        graphics::legend(legendpos,c("Estimate",colnames(ss)[se.pos]),pch=15,pt.cex=1.5,col=Col(c(density.col[1],se.col),alpha),lty=0,cex=0.6)
                    } else {
                        graphics::legend(legendpos,c("Estimate",colnames(ss)[se.pos]),col=c(density.col[1],se.col),lty=c(density.lty[1],se.lty),lwd=c(density.lwd[1],se.lwd),cex=0.6)
                }
                }
            }
            if (i==1 && ask) par(ask=ask)
        }
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
