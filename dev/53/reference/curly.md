# Adds curly brackets to plot

Adds curly brackets to plot

## Usage

``` r
curly(
  x,
  y,
  len = 1,
  theta = 0,
  wid,
  shape = 1,
  col = 1,
  lwd = 1,
  lty = 1,
  grid = FALSE,
  npoints = 50,
  text = NULL,
  offset = c(0.05, 0)
)
```

## Arguments

- x:

  center of the x axis of the curly brackets (or start end coordinates
  (x1,x2))

- y:

  center of the y axis of the curly brackets (or start end coordinates
  (y1,y2))

- len:

  Length of the curly brackets

- theta:

  angle (in radians) of the curly brackets orientation

- wid:

  Width of the curly brackets

- shape:

  shape (curvature)

- col:

  color (passed to lines/grid.lines)

- lwd:

  line width (passed to lines/grid.lines)

- lty:

  line type (passed to lines/grid.lines)

- grid:

  If TRUE use grid graphics (compatability with ggplot2)

- npoints:

  Number of points used in curves

- text:

  Label

- offset:

  Label offset (x,y)

## Examples

``` r
if (interactive()) {
plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
curly(x=c(1,0),y=c(0,1),lwd=2,text="a")
curly(x=c(1,0),y=c(0,1),lwd=2,text="b",theta=pi)
curly(x=-0.5,y=0,shape=1,theta=pi,text="c")
curly(x=0,y=0,shape=1,theta=0,text="d")
curly(x=0.5,y=0,len=0.2,theta=pi/2,col="blue",lty=2)
curly(x=0.5,y=-0.5,len=0.2,theta=-pi/2,col="red",shape=1e3,text="e")
}
```
