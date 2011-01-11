\name{compare}
\alias{compare}
\title{ Statistical tests }
\description{
  Performs Likelihood-ratio, Wald or score tests.
}
\usage{
compare(object, ...)
}
\arguments{
  \item{object}{ \code{lvmfit}-object }
  \item{\dots}{Additional arguments to low-level functions}
}
\value{
  Matrix of test-statistics and p-values
}
\examples{
m <- lvm(); 
regression(m) <- c(y1,y2,y3) ~ eta; latent(m) <- ~eta
regression(m) <- eta ~ x
}
\author{ Klaus K. Holst }
\seealso{ \code{\link{modelsearch}},  \code{\link{equivalence}} }
\keyword{ htest }
