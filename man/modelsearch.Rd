\name{modelsearch}
\alias{modelsearch}
\title{ Inference/model searching }
\description{
  Performs Wald or score tests.
}
\usage{
modelsearch(x, k = 1, dir = "forward")
}
\arguments{
  \item{x}{ \code{lvmfit}-object }
  \item{k}{ Number of parameters to test simultaneously }
  \item{dir}{ Direction to do model search. "forward" := add
    associations/arrows to model/graph (score tests), "backward" :=
    remove associations/arrows from model/graph (wald test) }
}
\value{
  Matrix of test-statistics and p-values
}
\examples{
m <- lvm(); 
regression(m) <- c(y1,y2,y3) ~ eta; latent(m) <- ~eta
regression(m) <- eta ~ x
m0 <- m; regression(m0) <- y2 ~ x
dd <- sim(m0,100)[,manifest(m0)]
e <- estimate(m,dd);
modelsearch(e)
}
\author{ Klaus K. Holst }
\seealso{ \code{\link{compare}}, }
\keyword{ htest }
