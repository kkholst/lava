##' Plot method for 'estimate' objects
##'
##' Plot method for 'estimate' objects
##' @export
##' @param x estimate object
##' @param f function of parameter coefficients and data parsed on to 'estimate'.
##'          If omitted a forest-plot will be produced.
##' @param idx Index of parameters (default all)
##' @param intercept include intercept in forest-plot
##' @param data data.frame
##' @param confint Add confidence limits
##' @param type plot type ('l')
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param col color
##' @param add add plot to current device
##' @param ... additional arguments to lower-level functions
plot.estimate <- function(x,f,idx,intercept=FALSE,data,confint=TRUE,type="l",xlab="x",ylab="f(x)",col=1,add=FALSE,...) {
    if (!missing(f) && !is.null(f)) {
        data <- as.list(data)
        env <- new.env()
        for (y in names(data)) {
            assign(y,data[[y]],env)
        }
        environment(f) <- env
        pp <- estimate(x,f,..., vcov=vcov(x),iid=FALSE)$coefmat
        if (!add) suppressWarnings(plot(data[[1]],pp[,1],xlab=xlab,ylab=ylab,type=type,...))
        else lines(data[[1]],pp[,1],xlab=xlab,ylab=ylab,type=type,col=col,...)
        if (confint) confband(data[[1]],pp[,3],pp[,4],polygon=TRUE,col=Col(col),lty=0)
        return(invisible(pp))
    }
    if (!is.null(x$coefmat)) {
        pp <- x$coefmat[,c(1,3,4),drop=FALSE]
    } else {
        pp <- cbind(coef(x),confint(x))
    }
    if (!missing(idx)) pp <- pp[idx,,drop=FALSE]
    if (!intercept) {
        idx <- match("(Intercept)",rownames(pp))
        if (length(idx)>0 && !any(is.na(idx))) pp <- pp[-idx,,drop=FALSE]
    }
    forestplot(pp[rev(seq(nrow(pp))),,drop=FALSE],...)
}

