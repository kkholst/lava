##' Adds curly brackets to plot
##'
##' @title Adds curly brackets to plot
##' @param x center of the x axis of the curly brackets (or start end coordinates (x1,x2))
##' @param y center of the y axis of the curly brackets (or start end coordinates (y1,y2))
##' @param len Length of the curly brackets
##' @param theta angle (in radians) of the curly brackets orientation
##' @param wid Width of the curly brackets
##' @param shape shape (curvature)
##' @param col color (passed to lines/grid.lines)
##' @param lwd line width (passed to lines/grid.lines)
##' @param lty line type (passed to lines/grid.lines)
##' @param grid If TRUE use grid graphics (compatability with ggplot2)
##' @param npoints Number of points used in curves
##' @param text Label
##' @param offset Label offset (x,y)
##' @export
##' @examples
##' if (interactive()) {
##' plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
##' curly(x=c(1,0),y=c(0,1),lwd=2,text="a")
##' curly(x=c(1,0),y=c(0,1),lwd=2,text="b",theta=pi)
##' curly(x=-0.5,y=0,shape=1,theta=pi,text="c")
##' curly(x=0,y=0,shape=1,theta=0,text="d")
##' curly(x=0.5,y=0,len=0.2,theta=pi/2,col="blue",lty=2)
##' curly(x=0.5,y=-0.5,len=0.2,theta=-pi/2,col="red",shape=1e3,text="e")
##' }
curly <- function(x,y,len=1,theta=0,
           wid,shape=1,
           col=1,lwd=1,lty=1,
           grid=FALSE,npoints=50,text=NULL,offset=c(0.05,0)) {
    if (length(x)==2 || length(y)==2) {
        x <- rep(x,length.out=2)
        y <- rep(y,length.out=2)
        v <- c(x[1]-x[2],y[1]-y[2])
        v0 <- c(1,0)-v
        len <- sum(v^2)^.5
        innerprod <- sum(v0)
        theta <- acos(innerprod/len)+theta
        len <- len/2
        x <- (x[1]-x[2])/2
        y <- (y[2]-y[1])/2
    }
    ii <- seq(0, pi/2, length.out=npoints)
    if (missing(wid)) {
        wid <- with(devcoords(),
        (fig.y2-fig.y1)/50)
    }
    x1 <- c(wid*(sin(ii)-1),
            c(0,0),
            wid*(1 - sin(rev(ii))),
            wid*(1 - sin(ii)),
            c(0,0),
            wid*(sin(rev(ii)) - 1))
    y1 <- c(-cos(ii),
            c(0,shape),
            shape+(cos(rev(ii))),
            shape+(2 - cos(ii)),
            c(shape+2, 2*shape+2),
            2*shape+2+cos(rev(ii)))

    x1 <- x1 + x + wid
    idx.max <- which.max(x1)
    max.y <- max(y1)
    y1 <- y1+1-(max.y+1)/2
    min.y <- min(y1)
    y1 <- y1*len/min.y+y
    ## Rotation
    x2 <- cos(theta) * (x1 - x) - sin(theta) * (y1 - y) + x
    y2 <- cos(theta) * (y1 - y) + sin(theta) * (x1 - x) + y
    x0 <- x1[idx.max]+offset[1]
    y0 <- y1[idx.max]+offset[2]
    xm <- cos(theta) * (x0 - x) - sin(theta) * (y0 - y) + x
    ym <- cos(theta) * (y0 - y) + sin(theta) * (x0 - x) + y
    if(grid){
        grid::grid.lines(grid::unit(x2,"npc"), grid::unit(y2,"npc"),
                         gp=grid::gpar(col=col,lwd=lwd,lty=lty))
    }
    else{
        points(x2,y2,type='l',col=col,lwd=lwd,lty=lty,xpd=TRUE)
    }
    theta <- acos(abs(cos(theta)))
    deg <- ((theta-pi/2)*180/pi)
    if (!is.null(text)) {
        text(xm,ym,text,srt=deg)
    }
}

