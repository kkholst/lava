\name{modelsearch}
\alias{modelsearch}
\alias{equivalence}
\title{ Model searching }
\description{
  Performs Wald or score tests and identifies candidates of equivalent models
}
\usage{
modelsearch(x, k = 1, dir = "forward",...)
equivalence(x, rel, tol = 0.001, k = 1, omitrel = TRUE, ...)
}
\arguments{
  \item{x}{ \code{lvmfit}-object }
  \item{k}{ Number of parameters to test simultaneously. For
    \code{equivalence} the number of additional associations to be added
  instead of \code{rel}.}
  \item{dir}{ Direction to do model search. "forward" := add
    associations/arrows to model/graph (score tests), "backward" :=
    remove associations/arrows from model/graph (wald test) }
  \item{rel}{Formula or character-vector specifying two variables to
    omit from the model and subsequently search for possible equivalent
    models}
  \item{tol}{ Define two models as empirical equivalent if the
    absolute difference in score test is less than \code{tol}}
  \item{omitrel}{if \code{k} greater than 1, this boolean defines wether
  to omit candidates containing \code{rel} from the output}
  \item{\dots}{Additional arguments to be passed to the low level functions}    
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
