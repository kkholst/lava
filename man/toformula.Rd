\name{toformula}
\alias{toformula}
\title{Converts strings to formula}
\usage{
  toformula(y = ".", x = ".")
}
\arguments{
  \item{y}{vector of predictors}

  \item{x}{vector of responses}
}
\value{
  An object of class \code{formula}
}
\description{
  Converts a vector of predictors and a vector of responses
  (characters) i#nto a formula expression.
}
\examples{
toformula(c("age","gender"), "weight")
}
\author{
  Klaus K. Holst
}
\seealso{
  \code{\link{as.formula}},
}
\keyword{models}
\keyword{utilities}

