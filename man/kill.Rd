\name{kill}
\alias{kill}
\alias{kill<-}
\title{Remove variables from (model) object.}
\usage{
  kill(x, ...) <- value
}
\arguments{
  \item{x}{Model object}

  \item{value}{Vector of variables or formula specifying
  which nodes to remove/add (or associations between
  nodes).}

  \item{\dots}{..}
}
\description{
  Generic method for removing elements of object
}
\examples{
m <- lvm()
addvar(m) <- ~y1+y2+x
covariance(m) <- y1~y2
regression(m) <- c(y1,y2) ~ x
### Cancel the covariance between the residuals of y1 and y2
cancel(m) <- y1~y2
### Remove y2 from the model
kill(m) <- ~y2
}
\author{
  Klaus K. Holst
}
\seealso{
  \code{cancel}
}
\keyword{models}
\keyword{regression}

