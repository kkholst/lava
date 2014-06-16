##' Extension of the \code{identify} function
##'
##' For the usual 'X11' device the identification process is
##' terminated by pressing any mouse button other than the first. For
##' the 'quartz' device the process is terminated by pressing either
##' the pop-up menu equivalent (usually second mouse button or
##' 'Ctrl'-click) or the 'ESC' key.
##' @title Identify points on plot
##' @usage 
##' \method{click}{default}(x, y=NULL, label=TRUE, n=length(x), pch=19, col="orange", cex=3, ...)
##' idplot(x,y,...,id=list())
##' @aliases idplot click.default
##' @param x X coordinates
##' @param y Y coordinates
##' @param label Should labels be added?
##' @param n Max number of inputs to expect
##' @param pch Symbol
##' @param col Colour
##' @param cex Size
##' @param id List of arguments parsed to \code{click} function
##' @param \dots Additional arguments parsed to \code{plot} function
##' @author Klaus K. Holst
##' @seealso \code{\link{idplot}}, \code{identify}
##' @examples
##' \donttest{
##' n <- 10; x <- seq(1:n); y <- runif(n)
##' plot(y ~ x); click(x,y)
##' 
##' data(iris)
##' l <- lm(Sepal.Length ~ Sepal.Width*Species,iris)
##' res <- plotConf(l,var2="Species")## ylim=c(6,8), xlim=c(2.5,3.3))
##' with(res, click(x,y))
##' 
##' with(iris, idplot(Sepal.Length,Petal.Length))
##' }
##' @keywords iplot
##' @export
click <- function(x,...){
    UseMethod("click")
}

##' @export
click.default <-
function(x, y=NULL, label=TRUE, n=length(x), pch=19, col="orange", cex=3, ...)
  {
    xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    sel <- rep(FALSE, length(x)); res <- integer(0)
    while(sum(sel) < n) {
      ans <- identify(x[!sel], y[!sel], n=1, plot=FALSE, ...)
      if(!length(ans)) break
      ans <- which(!sel)[ans]
             points(x[ans], y[ans], pch = pch, col=col, cex=cex)
      if (label)
        text(x[ans], y[ans], ans)
      sel[ans] <- TRUE
      res <- c(res, ans)
    }
    res

  }

##' @export
idplot <- function(x,y,...,id=list()) {
  plot(x,y,...)
  id$x <- x; id$y <- y
  do.call("click",id)
}
