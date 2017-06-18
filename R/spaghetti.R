##' Spaghetti plot for longitudinal data
##'
##' @title Spaghetti plot
##' @param formula Formula (response ~ time)
##' @param data data.frame
##' @param id Id variable
##' @param group group variable
##' @param type Type (line 'l', stair 's', ...)
##' @param lty Line type
##' @param pch Colour
##' @param col Colour
##' @param alpha transparency (0-1)
##' @param lwd Line width
##' @param level Confidence level
##' @param trend.formula Formula for trendline
##' @param tau Quantile to estimate (trend)
##' @param trend.lty Trend line type
##' @param trend.join Trend polygon
##' @param trend.delta Length of limit bars
##' @param trend Add trend line
##' @param trend.col Colour of trend line
##' @param trend.alpha Transparency
##' @param trend.lwd Trend line width
##' @param trend.jitter Jitter amount
##' @param legend Legend
##' @param by make separate plot for each level in 'by' (formula, name of column, or vector)
##' @param xlab Label of X-axis
##' @param ylab Label of Y-axis
##' @param add Add to existing device
##' @param ... Additional arguments to lower level arguments
##' @author Klaus K. Holst
##' @export
##' @examples
##' if (interactive() & requireNamespace("mets")) {
##' K <- 5
##' y <- "y"%++%seq(K)
##' m <- lvm()
##' regression(m,y=y,x=~u) <- 1
##' regression(m,y=y,x=~s) <- seq(K)-1
##' regression(m,y=y,x=~x) <- "b"
##' N <- 50
##' d <- sim(m,N); d$z <- rbinom(N,1,0.5)
##' dd <- mets::fast.reshape(d); dd$num <- dd$num+3
##' spaghetti(y~num,dd,id="id",lty=1,col=Col(1,.4),
##'           trend.formula=~factor(num),trend=TRUE,trend.col="darkblue")
##' dd$num <- dd$num+rnorm(nrow(dd),sd=0.5) ## Unbalance
##' spaghetti(y~num,dd,id="id",lty=1,col=Col(1,.4),
##'           trend=TRUE,trend.col="darkblue")
##' spaghetti(y~num,dd,id="id",lty=1,col=Col(1,.4),
##'            trend.formula=~num+I(num^2),trend=TRUE,trend.col="darkblue")
##' }
spaghetti <- function(formula,data,id="id",group=NULL,
              type="o",lty=1,pch=NA,col=1:10,alpha=0.3,lwd=1,
              level=0.95,
              trend.formula=formula,tau=NULL,
              trend.lty=1,trend.join=TRUE,trend.delta=0.2,
              trend=!is.null(tau),trend.col=col,
              trend.alpha=0.2,trend.lwd=3,
              trend.jitter=0,
              legend=NULL, by=NULL,
              xlab="Time",ylab="",add=FALSE,...) {
    ##spaghetti <- function(formula,data,id,type="l",lty=1,col=Col(1),trend=FALSE,trend.col="darkblue",trend.alpha=0.2,trend.lwd=3,xlab="Time",ylab="",...) {
    if (!lava.options()$cluster.index) stop("mets not available? Check 'lava.options()cluster.index'.")
    if (!is.null(by)) {
        if (is.character(by) && length(by==1)) {
            by <- data[,by]
        } else if (inherits(by,"formula")) {
            ##by <- model.matrix(update(by,~-1+.), model.frame(~.,data,na.action=na.pass))
            by <- model.frame(by,data,na.action=na.pass)
        }
        cl <- match.call(expand.dots=TRUE)
        cl$by <- NULL
        datasets <- split(data,by)
        res <- c()
        for (d in datasets) {
            cl$data <- d
            res <- c(res, eval(cl,parent.frame()))
        }
        return(invisible(res))
    }
    if (!is.null(group)) {
        if (is.character(group) && length(group==1)) {
            M <- data[,group]
        } else if (inherits(group,"formula")) {
            ##M <- model.matrix(update(group,~-1+.),data)
            M <- model.frame(group,data,na.action=na.pass)
        } else {
            M <- group
        }
        if (!add) plot(formula,data=data,xlab=xlab,ylab=ylab,...,type="n")
        dd <- split(data,M)
        K <- length(dd)
        if (length(type)<K)        type <- rep(type,K)
        if (length(col)<K)         col <- rep(col,K)
        if (length(pch)<K)         pch <- rep(pch,K)
        if (length(lty)<K)         lty <- rep(lty,K)
        if (length(lwd)<K)         lwd <- rep(lwd,K)
        if (length(alpha)<K)       alpha <- rep(alpha,K)
        if (length(trend)<K)       trend <- rep(trend,K)
        if (length(trend.col)<K)   trend.col <- rep(trend.col,K)
        if (length(trend.lty)<K)   trend.lty <- rep(trend.lty,K)
        if (length(trend.alpha)<K) trend.alpha <- rep(trend.alpha,K)
        if (length(trend.lwd)<K)   trend.lwd <- rep(trend.lwd,K)
        for (i in seq_len(K)) {
            spaghetti(formula,data=dd[[i]],id=id,type=type[i],
                     lty=lty[i],pch=pch[i],col=col[i],lwd=lwd[i],
                     alpha=alpha[i],
                     group=NULL,
                     trend=trend[i],tau=tau,
                     trend.col=trend.col[i],
                     trend.alpha=trend.alpha[i],
                     trend.lwd=trend.lwd[i],
                     trend.lty=trend.lty[i],
                     trend.delta=trend.delta,
                     trend.formula=trend.formula,
                     add=TRUE,...)
        }
        if (!is.null(legend)) {
            graphics::legend(legend,names(dd),lwd=lwd,col=col,lty=lty)
        }
        return(invisible(NULL))
    }

    if (inherits(id,"formula")) id <- all.vars(id)
    if (inherits(group,"formula")) group <- all.vars(group)
    if (is.character(id) && length(id)==1) Id <-
    y <- getoutcome(formula)
    x <- attributes(y)$x
    Idx <- function(vname,widenames) {
        idx <- which(unlist(lapply(widenames,function(x) length(grep(vname,substr(x,1,nchar(vname))))>0)))
        nn <- widenames[idx]
        ord <- order(char2num(unlist(lapply(nn,function(x) gsub(vname,"",x)))))
        idx[ord]
    }

    if (!requireNamespace("mets",quietly=TRUE)) stop("'mets' package required")

    if (length(x)==0) {
        data <- data[,c(id,y),drop=FALSE]
        wide <- mets::fast.reshape(data,id=id,varying=y,...)
        yidx <- Idx(y,names(wide))
        Y <- wide[,yidx,drop=FALSE]
        X <- NULL
        matplot(t(Y),type=type,lty=lty,pch=pch,lwd=lwd,col=Col(col[1],alpha[1]),xlab=xlab,ylab=ylab,...)
    } else {
        data <- data[,c(id,x,y),drop=FALSE]
        wide <- mets::fast.reshape(data[order(data[,id],data[,x]),],id=id,varying=c(y,x),...)
        yidx <- Idx(y,names(wide))
        xidx <- Idx(x,names(wide))
        Y <- wide[,yidx,drop=FALSE]
        X <- wide[,xidx,drop=FALSE]
        matplot(t(X),t(Y),type=type,pch=pch,lty=lty,lwd=lwd,col=Col(col[1],alpha[1]),xlab=xlab,ylab=ylab,add=add,...)
        if (trend) {
            if (is.numeric(trend.formula)) {
                trend.formula <- sort(trend.formula)
                tf <- toformula(y,"1")
                res <- c()
                if (!is.null(tau)) {
                    if (length(trend.alpha)<length(tau)) 	trend.alpha <- rep(trend.alpha,length(tau))
                    if (length(trend.lty)<length(tau)) 		trend.lty <- rep(trend.lty,length(tau))
                    if (length(trend.col)<length(tau)) 		trend.col <- rep(trend.col,length(tau))
                    if (length(trend.lwd)<length(tau)) 		trend.lwd <- rep(trend.lwd,length(tau))
                }
                for (i in trend.formula) {
                    data0 <- data[data[,x]==i,,drop=FALSE]
                    newdata <- data.frame(i); names(newdata) <- x
                    if (!is.null(tau)) {
                        ##if (!require(quantreg)) stop("Install 'quantreg'")
                        suppressWarnings(r1 <- quantreg::rq(tf,data=data0,tau=tau))
                        pr <- predict(r1,newdata=newdata,level=level)
                        res <- rbind(res,pr)
                    } else {
                        l1 <- lm(tf,data0)
                        pr <- predict(l1,newdata=newdata,interval="confidence",level=level)
                        res <- rbind(res,pr)
                    }
                }
                if (!is.null(tau)) {
                    for (j in seq_len(ncol(res))) {
                            if (trend.join) lines(trend.formula,res[,j],col=trend.col[j],lwd=trend.lwd[j],lty=trend.lty[j],...)
                            if (trend.delta>0) confband(trend.formula,res[,j],line=FALSE,col=trend.col[j],lty=trend.lty[j],lwd=trend.lwd[j],delta=trend.delta,...)
                        }
                } else {
                    confband(trend.formula,res[,2],res[,3],res[,1],col=Col(trend.col,trend.alpha),lty=trend.lty,lwd=trend.lwd,polygon=trend.join,...)
                }
            } else {
                tf <- getoutcome(trend.formula)
                if (is.list(tf)) {
                    trend.formula <- update(trend.formula,toformula(y,"."))
                }
                if (!is.null(tau)) {
                    ##if (!require(quantreg)) stop("Install 'quantreg'")
                    suppressWarnings(r1 <- quantreg::rq(trend.formula,data=data,tau=tau))
                    newdata <- data.frame(seq(min(X,na.rm=TRUE),max(X,na.rm=TRUE),length.out=100))
                    names(newdata) <- x
                    pr <- predict(r1,newdata=newdata,interval="confidence",level=level)
                    ##confband(xx,pr[,3],pr[,2],polygon=TRUE,col=Col(trend.col,trend.alpha),border=FALSE)
                    for (i in seq_along(tau))
                        lines(newdata[,1],pr[,i],col=trend.col,lwd=trend.lwd,lty=trend.lty)
                } else {
                    l1. <- lm(trend.formula,data)
                    l1 <- estimate(l1.,id=data[,id],level=level)
                    xy <- plotConf(l1.,vcov=vcov(l1),data=data,partres=FALSE,plot=FALSE,level=level,...)
                    xx <- xy$x
                    pr <- xy$predict$fit
                    if (is.factor(xx)) {
                        xx <- char2num(as.character(xx))
                        if (trend.jitter>0) xx <- jitter(xx,trend.jitter)                        
                        confband(xx,pr[,3],pr[,2],pr[,1],col=trend.col,lwd=2)
                    } else {
                        confband(xx,pr[,3],pr[,2],polygon=TRUE,col=Col(trend.col,trend.alpha),border=FALSE)
                        lines(xx,pr[,1],col=trend.col,lwd=trend.lwd,lty=trend.lty)
                    }
                }
            }
        }
    }
    return(invisible(list(Y,X)))
}

