##' Spaghetti plot for longitudinal data
##'
##' @title Spaghetti plot
##' @param formula Formula (response ~ time) 
##' @param data data.frame
##' @param id Id variable 
##' @param type Type (line 'l', stair 's', ...)
##' @param lty Line type
##' @param col Colour
##' @param trend Add trend line
##' @param trend.col Colour of trend line
##' @param trend.alpha Transparency
##' @param trend.lwd Trend line width
##' @param xlab Label of X-axis
##' @param ylab Label of Y-axis
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
##' d <- sim(m,500)
##' dd <- mets::fast.reshape(d);
##' dd$num <- dd$num+rnorm(nrow(dd),sd=0.5) ## Unbalance
##' spaghetti(y~num,dd,id="id",lty=1,col=Col(1,.4),trend=TRUE,trend.col="darkblue")
##' }
spaghetti <- function(formula,data,id,type="l",lty=1,col=Col(1),trend=FALSE,trend.col="darkblue",trend.alpha=0.2,trend.lwd=3,xlab="Time",ylab="",...) {
    if (!lava.options()$cluster.index) stop("mets not available? Check 'lava.options()cluster.index'.")
    y <- getoutcome(formula)
    x <- attributes(y)$x
    Idx <- function(vname,widenames) {
        idx <- which(unlist(lapply(widenames,function(x) length(grep(vname,substr(x,1,nchar(vname))))>0)))
        nn <- widenames[idx]
        ord <- order(as.numeric(unlist(lapply(nn,function(x) gsub(vname,"",x)))))
        idx[ord]
    }
    if (length(col)==nrow(data)) {
        col <- with(data, fast.reshape(col,id=id))[,1]
    }
    if (length(x)==0) {
        wide <- mets::fast.reshape(data,id=id,varying=y,...)
        yidx <- Idx(y,names(wide))
        Y <- wide[,yidx,drop=FALSE]
        X <- NULL
        matplot(t(Y),type=type,lty=lty,col=col,xlab=xlab,ylab=ylab,...)
    } else {
        wide <- mets::fast.reshape(dsort(data,id,x),id=id,varying=c(y,x),...)
        yidx <- Idx(y,names(wide))
        xidx <- Idx(x,names(wide))
        Y <- wide[,yidx,drop=FALSE]
        X <- wide[,xidx,drop=FALSE]
        matplot(t(X),t(Y),type=type,lty=lty,col=col,xlab=xlab,ylab=ylab,...)
        if (trend) {
            l1. <- lm(formula,data)
            l1 <- estimate(l1.,id=data[,id])
            xy <- plotConf(l1.,vcov=vcov(l1),partres=FALSE,plot=FALSE,...)
            xx <- xy$predict.newdata[,x]
            pr <- xy$predict$fit
            confband(xx,pr[,3],pr[,2],polygon=TRUE,col=Col(trend.col,trend.alpha))
            lines(xx,pr[,1],col=trend.col,lwd=trend.lwd)
        }
    }
    return(invisible(list(Y,X)))
}
