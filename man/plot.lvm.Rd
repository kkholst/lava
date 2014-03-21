\name{plot.lvm}
\alias{plot.lvm}
\alias{plot.lvmfit}
\title{Plot path diagram}
\usage{
  \method{plot}{lvm} (x, diag = FALSE, cor = TRUE,
    labels = FALSE, intercept = FALSE, addcolor = TRUE,
    plain = FALSE, cex, fontsize1 = 10, noplot = FALSE,
    graph = list(rankdir = "BT"),
    attrs = list(graph = graph), unexpr = FALSE,
    addstyle = TRUE, Rgraphviz = lava.options()$Rgraphviz,
    init = TRUE,
    layout = c("dot", "fdp", "circo", "twopi", "neato", "osage"),
    ...)
}
\arguments{
  \item{x}{Model object}

  \item{diag}{Logical argument indicating whether to
  visualize variance parameters (i.e. diagonal of variance
  matrix)}

  \item{cor}{Logical argument indicating whether to
  visualize correlation parameters}

  \item{labels}{Logical argument indiciating whether to add
  labels to plot (Unnamed parameters will be labeled
  p1,p2,...)}

  \item{intercept}{Logical argument indiciating whether to
  add intercept labels (current version: not used))}

  \item{addcolor}{Logical argument indiciating whether to
  add colors to plot (overrides \code{nodecolor} calls)}

  \item{plain}{if TRUE strip plot of colors and boxes}

  \item{cex}{Fontsize of node labels}

  \item{fontsize1}{Fontsize of edge labels}

  \item{noplot}{if TRUE then return \code{graphNEL} object
  only}

  \item{graph}{Graph attributes (Rgraphviz)}

  \item{attrs}{Attributes (Rgraphviz)}

  \item{unexpr}{if TRUE remove expressions from labels}

  \item{addstyle}{Logical argument indicating whether
  additional style should automatically be added to the
  plot (e.g. dashed lines to double-headed arrows)}

  \item{Rgraphviz}{if FALSE igraph is used for graphics}

  \item{init}{Reinitialize graph (for internal use)}

  \item{layout}{Graph layout (see Rgraphviz or igraph
  manual)}

  \item{\dots}{Additional arguments to be passed to the low
  level functions}
}
\description{
  Plot the path diagram of a SEM
}
\examples{
\donttest{
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
}

\donttest{
m <- lvm(c(y1,y2) ~ eta)
regression(m) <- eta ~ z+x2
regression(m) <- c(eta,z) ~ x1
latent(m) <- ~eta
plot(lava:::beautify(m,edgecol=F))
}
}
\author{
  Klaus K. Holst
}
\keyword{hplot}
\keyword{regression}

