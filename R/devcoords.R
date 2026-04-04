##' Returns device-coordinates and plot-region
##'
##' @title Returns device-coordinates and plot-region
##' @return A `list` with elements
##'  * `dev.x1`: device left x-coordinate
##'  * `dev.x2`: device right x-coordinate
##'  * `dev.y1`: device bottom y-coordinate
##'  * `dev.y2`: device top y-coordinate
##'  * `fig.x1`: plot left x-coordinate
##'  * `fig.x2`: plot right x-coordinate
##'  * `fig.y1`: plot bottom y-coordinate
##'  * `fig.y2`: plot top y-coordinate
##' @author Klaus K. Holst
##' @export
##' @keywords hplot
`devcoords` <- function() {
  cc <- par("usr") ## extremes of coordinates of plotting region (x1,x2,y1,y2)
  plotinch <- par("pin") ## Plot dimensions (width,height) in inches
  margininch <- par("mai") ## Margin sizes in inches (bottom, left, top ,right)
  plotlenX <- cc[2]-cc[1]
  unitinchX <- plotlenX/plotinch[1]
  plotlenY <- cc[4]-cc[3]
  unitinchY <- plotlenY/plotinch[2]
  deviceXleft <- cc[1]-unitinchX*margininch[2]
  deviceXright <- cc[2]+unitinchX*margininch[4]
  deviceYtop <- cc[4]+unitinchY*margininch[3]
  deviceYbottom <- cc[3]-unitinchY*margininch[1]
  return(list(dev.x1=deviceXleft, dev.x2=deviceXright,
              dev.y1=deviceYbottom, dev.y2=deviceYtop,
              fig.x1=cc[1], fig.x2=cc[2],
              fig.y1=cc[3], fig.y2=cc[4]))
}
