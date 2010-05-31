\name{plot.lvm}
\alias{plot.lvm}
\alias{plot.lvmfit}
\title{ Plot path diagram }
\description{
  Plot the path diagram of a SEM
}
\usage{
\method{plot}{lvm} (x, diag = FALSE, cor=TRUE, labels = FALSE,
addcolor = TRUE, debug=FALSE, ...)
}
\arguments{
  \item{x}{ Model object }
  \item{diag}{ Logical argument indicating whether to visualize variance
    parameters (i.e. diagonal of variance matrix)}
  \item{cor}{ Logical argument indicating whether to visualize
    correlation parameters }
  \item{labels}{ Logical argument indiciating whether to add labels to
    plot (Unnamed parameters will be labeled p1,p2,...)}
  \item{addcolor}{ Logical argument indiciating whether to add colors to
    plot (overrides \code{nodecolor} calls) }
  \item{debug}{ for debugging only }
  \item{\dots}{}
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
