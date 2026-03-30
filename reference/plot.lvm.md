# Plot path diagram

Plot the path diagram of a SEM

## Usage

``` r
# S3 method for class 'lvm'
plot(
  x,
  diag = FALSE,
  cor = TRUE,
  labels = FALSE,
  intercept = FALSE,
  addcolor = TRUE,
  plain = FALSE,
  cex,
  fontsize1 = 10,
  noplot = FALSE,
  graph = list(rankdir = "BT"),
  attrs = list(graph = graph),
  unexpr = FALSE,
  addstyle = TRUE,
  plot.engine = lava.options()$plot.engine,
  init = TRUE,
  layout = lava.options()$layout,
  edgecolor = lava.options()$edgecolor,
  graph.proc = lava.options()$graph.proc,
  ...
)
```

## Arguments

- x:

  Model object

- diag:

  Logical argument indicating whether to visualize variance parameters
  (i.e. diagonal of variance matrix)

- cor:

  Logical argument indicating whether to visualize correlation
  parameters

- labels:

  Logical argument indiciating whether to add labels to plot (Unnamed
  parameters will be labeled p1,p2,...)

- intercept:

  Logical argument indiciating whether to add intercept labels

- addcolor:

  Logical argument indiciating whether to add colors to plot (overrides
  `nodecolor` calls)

- plain:

  if TRUE strip plot of colors and boxes

- cex:

  Fontsize of node labels

- fontsize1:

  Fontsize of edge labels

- noplot:

  if TRUE then return `graphNEL` object only

- graph:

  Graph attributes (Rgraphviz)

- attrs:

  Attributes (Rgraphviz)

- unexpr:

  if TRUE remove expressions from labels

- addstyle:

  Logical argument indicating whether additional style should
  automatically be added to the plot (e.g. dashed lines to double-headed
  arrows)

- plot.engine:

  default 'Rgraphviz' if available, otherwise visNetwork,igraph

- init:

  Reinitialize graph (for internal use)

- layout:

  Graph layout (see Rgraphviz or igraph manual)

- edgecolor:

  if TRUE plot style with colored edges

- graph.proc:

  Function that post-process the graph object (default: subscripts are
  automatically added to labels of the nodes)

- ...:

  Additional arguments to be passed to the low level functions

## Author

Klaus K. Holst

## Examples

``` r
if (interactive()) {
m <- lvm(c(y1,y2) ~ eta)
regression(m) <- eta ~ z+x2
regression(m) <- c(eta,z) ~ x1
latent(m) <- ~eta
labels(m) <- c(y1=expression(y[scriptscriptstyle(1)]),
y2=expression(y[scriptscriptstyle(2)]),
x1=expression(x[scriptscriptstyle(1)]),
x2=expression(x[scriptscriptstyle(2)]),
eta=expression(eta))
edgelabels(m, eta ~ z+x1+x2, cex=2, lwd=3,
           col=c("orange","lightblue","lightblue")) <- expression(rho,phi,psi)
nodecolor(m, vars(m), border="white", labcol="darkblue") <- NA
nodecolor(m, ~y1+y2+z, labcol=c("white","white","black")) <- NA
plot(m,cex=1.5)

d <- sim(m,100)
e <- estimate(m,d)
plot(e)

m <- lvm(c(y1,y2) ~ eta)
regression(m) <- eta ~ z+x2
regression(m) <- c(eta,z) ~ x1
latent(m) <- ~eta
plot(lava:::beautify(m,edgecol=FALSE))
}
```
