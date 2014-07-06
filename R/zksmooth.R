##' Plot/estimate surface 
##'
##' @export
##' @aliases ksmooth surface
##' @param x formula or data
##' @param data data.frame
##' @param h bandwidth
##' @param rgl If TRUE 'rgl' is used
##' @param xlab X label
##' @param ylab Y label
##' @param zlab Z label
##' @param gridsize grid size of kernel smoother
##' @param ... Additional arguments to graphics routine (persp3d or persp)
##' @examples
##' 
##' f <- function(x,y) 1-sqrt(x^2+y^2)
##' surface(f,xlim=c(-1,1),alpha=0.9,aspect=c(1,1,0.75))
##' surface(f,xlim=c(-1,1),clut=heat.colors(128))
##' ##play3d(spin3d(axis=c(0,0,1), rpm=8), duration=5)
##' 
##' 
##' surface(function(x) dmvn(x,sigma=diag(2)),c(-3,3),lit=FALSE,smooth=FALSE,box=FALSE,alpha=0.8)
ksmooth <- function(x,data,h,rgl=TRUE,xlab,ylab,zlab="",gridsize=rep(51L,2),...) {
    if (inherits(x,"formula")) {
        x <- model.frame(x,data)
    }
    if (missing(h)) h <- sd(as.matrix(x))*nrow(x)^(-1/5)
    est <- KernSmooth::bkde2D(x, bandwidth=h, gridsize=gridsize)
    if (missing(xlab)) xlab <- names(x)[1]
    if (missing(ylab)) ylab <- names(x)[2]
    if (rgl) {
        surface(est$fhat, x=est$x1, y=est$x2, est$fhat, 
                xlab=xlab, ylab=ylab, zlab=zlab, ...)
    } else {
        op <- par(mfrow=c(2,1))
        persp(est$fhat, expand = 0.5, col=col, ...) ##theta = 15, phi = 25, 
        contour(est$x1, est$x2, est$fhat, ...) ##nlevels=20
        par(op)
    }
    return(invisible(est))
}


##' @export
surface <- function(f,xlim=c(0,1),ylim=xlim,n=rep(100,2),col,clut,x,y,...) {
    if (missing(x)) {
        x <- seq(xlim[1],xlim[2],length.out=n[1])    
        y <- seq(ylim[1],ylim[2],length.out=n[2])        
    }
    if (is.function(f)) {
        xy <- as.matrix(expand.grid(x,y))
        
        if (inherits(try(f(c(x[1],y[1])),silent=TRUE),"try-error"))
            f <- f(xy[,1],xy[,2])
        else
            f <- f(xy)
    }
    zrg <- range(f) 
    zlen <- diff(zrg)
    if (missing(clut)) {
        jet.colors <- colorRampPalette( c("blue","red"),bias=1)
        clut <- jet.colors(128)
    }
    if (missing(col)) {
        col <- clut[round((length(clut)-1)*(f-zrg[1])/zlen)+1]
    }
    rgl::persp3d(x, y, f, col=col,...)
}
