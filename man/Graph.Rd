\name{Graph}
\alias{Graph}
\alias{Graph.lvmfit}
\alias{Graph.lvm}
\alias{Graph<-}
\alias{Graph<-.lvmfit}
\alias{Graph<-.lvm}


\title{ Extract graph }
\description{
  Extract or replace graph object
}
\usage{
Graph(x, ...)
Graph(x, ...) <- value
}
\arguments{
  \item{x}{ Model object }
  \item{value}{ New \code{graphNEL} object }
  \item{\dots}{}
}
\value{
  Returns a graph object (\code{graphNEL})
}
\author{Klaus K. Holst}
\seealso{ \code{\link{Model}} }
\examples{
m <- lvm(y~x)
Graph(m)
}

\keyword{ graphs }
\keyword{ models }
