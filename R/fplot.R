##' Faster plot via RGL
##' @title fplot
##' @export
##' @examples
##' if (interactive()) {
##' data(iris)
##' fplot(Sepal.Length ~ Petal.Length+Species, data=iris, size=2, type="s")
##' }
##' @param x X variable
##' @param y Y variable
##' @param z Z variable (optional)
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param ... additional arggument to lower-level plot functions
##' @param z.col color (use argument alpha to set transparency)
##' @param data data.frame
##' @param add if TRUE use current active device
##' @param aspect aspect ratio
##' @param zoom zoom level
fplot <- function(x,y,z=NULL,xlab,ylab,...,z.col=topo.colors(64), 
                  data=parent.frame(),add=FALSE,aspect=c(1,1),zoom=0.8) {
    if (!requireNamespace("rgl",quietly=TRUE)) stop("Requires 'rgl'")    
    if (inherits(x,"formula")) {
        y <- getoutcome(x)
        x <- attributes(y)$x
        if (length(x)>1) {
            z <- as.numeric(with(data, get(x[2])))
        }
        if (length(x)==0) {
            x <- seq(nrow(data))
            if (missing(xlab)) xlab <- "Index"
        } else {
            if (missing(xlab)) xlab <- x[1]
            x <- with(data, get(x[1]))
        }
        if (missing(ylab)) ylab <- y
        y <- with(data, get(y))
    } else {
        if (missing(y)) {
            y <- x
            if (missing(ylab)) ylab <- deparse(substitute(x))
            x <- seq(nrow(data))            
            if (missing(xlab)) xlab <- "Index"       
        } else {
            if (missing(xlab)) xlab <- deparse(substitute(x))
            if (missing(ylab)) ylab <- deparse(substitute(y))
        }
    }
    rgl::.check3d()
    if (!is.null(z)) {
        ncol <- length(z.col);
        glut <- approxfun(seq(min(z),max(z),length.out=ncol),seq(ncol))
        rgl::plot3d(x,y,0,col=z.col[round(glut(z))],xlab=xlab,ylab=ylab,add=add,...)
    } else {
        rgl::plot3d(x,y,0,xlab=xlab,ylab=ylab,add=add,...)
        
    }
    rgl::view3d(0,0,fov=0,zoom=zoom)
    rgl::aspect3d(c(aspect,1))
}


