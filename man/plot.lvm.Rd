\name{plot.lvm}
\alias{plot.lvm}
\alias{plot.lvmfit}
\title{ Plot path diagram }
\description{
  Plot the path diagram of a SEM
}
\usage{
\method{plot}{lvm}(x, diag = FALSE, cor = TRUE, labels =
                 FALSE, intercept = FALSE, addcolor = TRUE, plain =
                 FALSE, cex, fontsize1 = 10, noplot =
                 FALSE, graph = list(rankdir = "BT"), attrs =
                 list(graph = graph), unexpr = FALSE,
                 addstyle = TRUE, Rgraphviz=lava.options()$Rgraphviz, ...)
}
\arguments{
  \item{x}{ Model object }
  \item{diag}{ Logical argument indicating whether to visualize variance
    parameters (i.e. diagonal of variance matrix)}
  \item{cor}{ Logical argument indicating whether to visualize
    correlation parameters }
  \item{labels}{ Logical argument indiciating whether to add labels to
    plot (Unnamed parameters will be labeled p1,p2,...)}
  \item{intercept}{ Logical argument indiciating whether to add
    intercept labels (current version: not used))}
  \item{addcolor}{ Logical argument indiciating whether to add colors to
    plot (overrides \code{nodecolor} calls) }
  \item{plain}{ if TRUE strip plot of colors and boxes}
  \item{cex}{ Fontsize of node labels }
  \item{fontsize1}{ Fontsize of edge labels}
  \item{noplot}{if TRUE then return \code{graphNEL} object only}
  \item{graph}{Graph attributes (Rgraphviz)}
  \item{attrs}{Attributes (Rgraphviz)}
  \item{unexpr}{if TRUE remove expressions from labels}
  \item{addstyle}{Logical argument indicating whether additional style
    should automatically be added to the plot (e.g. dashed lines to
    double-headed arrows) }
  \item{Rgraphviz}{if FALSE igraph is used for graphics }
  \item{\dots}{Additional arguments to be passed to the low level functions}
}
\author{ Klaus K. Holst }
\examples{
\dontrun{
example(estimate)
plot(e)
}
}
\keyword{ hplot }
\keyword{ regression }
