# Define labels of graph

Alters labels of nodes and edges in the graph of a latent variable model

## Usage

``` r
# Default S3 method
labels(object, ...) <- value
# S3 method for class 'lvm'
edgelabels(object, to, ...) <- value
# Default S3 method
nodecolor(object, var=vars(object),
border, labcol, shape, lwd, ...) <- value
```

## Arguments

- object:

  `lvm`-object.

- ...:

  Additional arguments (`lwd`, `cex`, `col`, `labcol`), `border`.

- value:

  node label/edge label/color

- to:

  Formula specifying outcomes and predictors defining relevant edges.

- var:

  Formula or character vector specifying the nodes/variables to alter.

- border:

  Colors of borders

- labcol:

  Text label colors

- shape:

  Shape of node

- lwd:

  Line width of border

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm(c(y,v)~x+z)
regression(m) <- c(v,x)~z
labels(m) <- c(y=expression(psi), z=expression(zeta))
nodecolor(m,~y+z+x,border=c("white","white","black"),
          labcol="white", lwd=c(1,1,5),
          lty=c(1,2)) <-  c("orange","indianred","lightgreen")
edgelabels(m,y~z+x, cex=c(2,1.5), col=c("orange","black"),labcol="darkblue",
           arrowhead=c("tee","dot"),
           lwd=c(3,1)) <- expression(phi,rho)
edgelabels(m,c(v,x)~z, labcol="red", cex=0.8,arrowhead="none") <- 2
if (interactive()) {
    plot(m,addstyle=FALSE)
}

m <- lvm(y~x)
labels(m) <- list(x="multiple\nlines")
if (interactive()) {
op <- par(mfrow=c(1,2))
plot(m,plain=TRUE)
plot(m)
par(op)

d <- sim(m,100)
e <- estimate(m,d)
plot(e,type="sd")
}
```
