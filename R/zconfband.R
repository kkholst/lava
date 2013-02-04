##' Add Confidence limits bar to plot
##'
##' @title Add Confidence limits bar to plot
##' @param x Position (x-coordinate if vert=TRUE, y-coordinate otherwise)
##' @param lower Lower limit (if NULL no limits is added, and only the
##' center is drawn (if not NULL))
##' @param upper Upper limit
##' @param center Center point
##' @param delta Length of limit bars
##' @param centermark Length of center bar
##' @param pch Center symbol (if missing a line is drawn)
##' @param blank If TRUE a white ball is plotted before the center is
##' added to the plot
##' @param vert If TRUE a vertical bar is plotted. Otherwise a horizontal
##' bar is used
##' @param \dots Additional low level arguments (e.g. col, lwd, lty,...)
##' @seealso \code{confband}
##' @export
##' @keywords iplot
##' @examples
##' plot(0,0,type="n",xlab="",ylab="")
##' confband(0.5,-0.5,0.5,0,col="darkblue")
##' confband(0.8,-0.5,0.5,0,col="darkred",vert=FALSE,pch=1,cex=1.5)
##' @author Klaus K. Holst
confband <- function(x,lower,upper,center=NULL,delta=0.07,centermark=0.03,
                     pch,blank=TRUE,vert=TRUE,...) {
  if (vert) {
    if (!is.null(lower)) {
      segments(x,lower,x,upper,...)
      segments(x-delta,lower,x+delta,lower,...)
      segments(x-delta,upper,x+delta,upper,...)
    }
    if (!is.null(center)) {
      if (!missing(pch)) {
        if (blank)
          points(x,center,pch=16,col="white")
        points(x,center,pch=pch,...)
      } else {
        segments(x-centermark,center,x+centermark,center,...)
      }
    }
  } else {
  
      if (!is.null(lower)) {
        segments(lower,x,upper,x,...)
        segments(lower,x-delta,lower,x+delta,...)
        segments(upper,x-delta,upper,x+delta,...)
      }
      if (!is.null(center)) {
        if (!missing(pch)) {
          if (blank)
        points(center,x,pch=16,col="white")
          points(center,x,pch=pch,...)
        } else {
          segments(center,x-centermark,center,x+centermark,...)
        }
      }
    }
    invisible(c(x,lower,upper,center))
}
