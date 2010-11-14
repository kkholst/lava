\name{Model}
\alias{Model}
\alias{Model<-}
\alias{Model.lvmfit}
\alias{Model<-.lvmfit}
\alias{Model.multigroupfit}
\alias{Model<-.multigroupfit}

\title{ Extract model }
\description{
  Extract or replace model object
}
\usage{
Model(x, ...)
Model(x, ...) <- value
}
\arguments{
  \item{x}{ Fitted model }
  \item{value}{ New model object (e.g. \code{lvm} or \code{multigroup})}
  \item{\dots}{Additional arguments to be passed to the low level functions}
}
\value{
  Returns a model object (e.g. \code{lvm} or \code{multigroup})
}
\author{Klaus K. Holst}
\seealso{ \code{\link{Graph}} }
\examples{
m <- lvm(y~x)
e <- estimate(m, sim(m,100))
Model(e)
}

\keyword{ models }
